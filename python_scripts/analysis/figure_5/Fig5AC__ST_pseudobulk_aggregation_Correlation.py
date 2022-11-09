#!/usr/bin/env python
"""Calculate Pearson Correlation between Cytokines and their corresponding responder genes UMI-counts
    File name: Fig4A__ST_pseudobulk_aggregation_Correlation.py
    Author: Christina Hillig
    Date created: 23/11/2020
    Date last modified: 4/29/2021
    Python Version: 3.7
"""

from python_scripts.spatial_correlation import corr_statistics as corr_stats
from python_scripts.utils import gene_lists
from python_scripts.spatial_correlation.plots import plot_colorbar_legend
from python_scripts.spatial_correlation import helper_functions as ht

import scanpy as sc
import scipy.stats as stats
import numpy as np
import pandas as pd
import itertools
import os
from datetime import date

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  # for legend handle


fig_size = (6, 6)
xy_fontsize = 16
xy_ticks = 12
title_fontsize = 18
legend_fontsize = 14
text_fontsize = 14
fileformat = '.pdf'
size_multiplier = 6


def get_bulk_cytoresp_counts(adata, genes_dict, tissue_types=None):
    """Get per sample the counts of a cytokine and its responder genes

    Parameters
    ----------
    adata : annData
    genes_dict : dict
        keys are the cytokines and entries are the responder genes
    tissue_types : str, list

    Returns
    -------

    """

    if tissue_types:
        bool_col = adata.obs[tissue_types] == 1
        merged = bool_col.sum(axis=1)
        adata = adata[merged == 1]

    specimen_name = np.unique(adata.obs['specimen'])
    # get frequency
    col_names = list(genes_dict.keys())
    for gen in genes_dict.keys():
        col_names.append("_".join([gen, 'responder']))
    col_names.extend(['disease', 'n_IFNG', 'n_IL13', 'n_IL17A'])

    df_corr = pd.DataFrame(index=specimen_name, columns=col_names)

    for ind, specimen in enumerate(specimen_name):
        ad = adata[adata.obs['specimen'] == specimen]
        df_corr['disease'][ind] = np.unique(ad.obs['DISEASE'].values)[0]
        for cytokine_gene in genes_dict:
            # 1. Get counts
            if "counts" in adata.layers.keys():
                counts_gene = ad.layers['counts'][:, np.where(ad.var.index == cytokine_gene)[0]]
            else:
                counts_gene = ad.X[:, np.where(ad.var.index == cytokine_gene)[0]]

            # sum counts of cytokine found in sample xy and add to dataframe
            counts_cyto_summed = np.asarray(counts_gene[:, 0], dtype=int).sum(axis=0)
            n_cytospots = np.count_nonzero(np.asarray(counts_gene[:, 0], dtype=int))

            # Get counts of responder genes and sum them up
            counts_responders = 0
            no_resps_spot = 0
            for responder_gene in genes_dict[cytokine_gene]:
                if responder_gene in ad.var.index:
                    # 1. Get counts
                    if "counts" in ad.layers.keys():
                        counts_gene = ad.layers['counts'][:, np.where(ad.var.index == responder_gene)[0]]
                    else:
                        counts_gene = ad.X[:, np.where(ad.var.index == responder_gene)[0]]

                    # Sum up counts of responder genes for each cytokine and add to dataframe
                    counts_responders += np.asarray(counts_gene[:, 0], dtype=int).sum(axis=0)
                    no_resps_spot += len(np.where(counts_gene[:, 0] != 0)[0])

            # Add counts to dataframe
            df_corr.at[specimen, 'n_{}'.format(cytokine_gene)] = n_cytospots
            df_corr.at[specimen, cytokine_gene] = counts_cyto_summed
            df_corr.at[specimen, "_".join([cytokine_gene, 'responder'])] = counts_responders

    return df_corr


def weighted_transcripts_correlation_plot(df_counts_cytoresps, genes_dict, corr_method, save_folder):
    """
    Plot correlation of counts distribution for each cytokine over all samples

    :param df_counts_cytoresps:
    :param genes_dict:
    :param corr_method:
    :param save_folder:
    :return:
    """

    # df_counts_cytoresps = df_counts_cytoresps.replace({'0': np.nan, 0: np.nan})
    df_counts_cytoresps = df_counts_cytoresps.replace({np.nan: 0})

    cyto_combs = list(set(itertools.combinations(genes_dict.keys(), 2)))
    colors = {'LP': 'orange', 'AD': 'red', 'Pso': 'blue'}

    # Correlation between cytokine counts
    for cyto_combination in cyto_combs:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.grid(False)
        corr, pval = stats.pearsonr(df_counts_cytoresps[cyto_combination[0]], df_counts_cytoresps[cyto_combination[1]])

        ax.scatter(df_counts_cytoresps[cyto_combination[0]], df_counts_cytoresps[cyto_combination[1]],
                   c=df_counts_cytoresps['disease'].map(colors))
        ax.set_ylabel(cyto_combination[1])
        ax.set_xlabel(cyto_combination[0])

        # Add text: Correlation value and p-value
        ax.text(
            df_counts_cytoresps[cyto_combination[0]].max() / 2 - df_counts_cytoresps[cyto_combination[0]].max() / 10,
            df_counts_cytoresps[cyto_combination[1]].max(),
            'r = {:.2f}; p = {:.2e}'.format(corr, pval), fontstyle='italic', fontsize=text_fontsize)

        # add a legend
        handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=8) for k, v in
                   colors.items()]
        ax.legend(title='Disease', handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # Save figure
        fig.savefig(os.path.join(
            save_folder, 'Fig4A-C__{}_{}_bulkapproach{}'.format(cyto_combination[0], cyto_combination[1], fileformat)))
        plt.close()

    # Correlation between cytokine responder counts
    cytoresps_combs = list(set(itertools.combinations(['IFNG_responder', 'IL13_responder', 'IL17A_responder'], 2)))
    colors = {'LP': 'tab:blue', 'AD': 'tab:orange', 'Pso': 'tab:red'}
    text_pos = {'Pso': df_counts_cytoresps.max()[1],
                'AD': df_counts_cytoresps.max()[1] / 2,
                'LP': df_counts_cytoresps.min()[1]}
    for cyto_combination in cytoresps_combs:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.grid(False)
        for diag in np.unique(df_counts_cytoresps['disease']):
            df_counts_cytoresps_temp = df_counts_cytoresps[df_counts_cytoresps['disease'] == diag]
            corr, pval = stats.pearsonr(df_counts_cytoresps_temp[cyto_combination[0]],
                                        df_counts_cytoresps_temp[cyto_combination[1]])

            ax.scatter(df_counts_cytoresps_temp[cyto_combination[0]], df_counts_cytoresps_temp[cyto_combination[1]],
                       c=colors[diag])
            ax.set_ylabel(cyto_combination[1])
            ax.set_xlabel(cyto_combination[0])

            # Add text: Correlation value and p-value
            ax.text(text_pos[diag], df_counts_cytoresps_temp.max()[0] + 1,
                    r'{}: r = {:.2f}; p = {:.2e}'.format(diag, corr, pval),
                    fontstyle='italic', fontsize=10)

        # add a legend
        handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=8) for k, v in
                   colors.items()]
        ax.legend(title='Disease', handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # Save figure
        fig.savefig(os.path.join(
            save_folder, 'Fig4A-C__{}_{}_bulkapproach{}'.format(cyto_combination[0], cyto_combination[1], fileformat)))
        plt.close()

    # Save dataframe
    for cyto in genes_dict.keys():
        weighted_cytoname = "_".join(['weighted', cyto])
        df_counts_cytoresps[weighted_cytoname] = df_counts_cytoresps[cyto]
    df_counts_cytoresps.to_csv(os.path.join(save_folder, 'Bulk_.csv'))

    dict_corr = dict({'pearson': [], 'spearman': []})
    for cyto in genes_dict.keys():
        weighted_cytoname = "_".join(['weighted', cyto])

        for cyto_reps in genes_dict.keys():
            resp_name = "_".join([cyto_reps, 'responder'])
            # temp_df = df_counts_cytoresps[[cyto, resp_name, 'disease', 'marker']]
            temp_df = df_counts_cytoresps[[cyto, resp_name, 'disease', 'n_{}'.format(cyto)]]

            # Apply weights to cytokine counts and add it to df: Multiply by cyto transcripts
            weighted_cytocounts = temp_df[cyto]
            temp_df[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)

            # 2.1 Calculate weighted Pearson Correlation per disease
            dict_corr = corr_stats.get_correlation_stats(
                df=temp_df, resp_name=resp_name, cyto=cyto, weight_name=cyto, dict_corr=dict_corr)

            # Read out correlation value and p-value for each disease
            corr_pval = [a_tuple for a_tuple in dict_corr[corr_method]
                         if (a_tuple[0] == cyto) and (a_tuple[1] == resp_name)]

            # Save df in a new df and sort it - otherwise fill_between wont work
            xy_vals = temp_df[[resp_name, cyto]].copy()
            xy_vals = xy_vals.sort_values(by=resp_name)

            # Prepare data for linear regression fit
            x = np.asarray([xy_vals[resp_name].values]).T
            y = np.asarray([xy_vals[cyto].values]).T

            est_woweights_fitted, ypred_woweights, df_woweights = corr_stats.create_ls_model(x=x, y=y, w=None)

            # Plot Correlation
            # Increase point size
            xy_vals["Cluster size"] = (xy_vals[cyto].values ** 2 * size_multiplier + 1).astype(int)
            #  2. Plot correlation
            fig, ax = plt.subplots(figsize=fig_size)
            ax.grid(False)
            # test = regr_wweight.coef_[0] * x + regr_wweight.intercept_
            # ax.plot(x, test, color='green', label='test', linewidth=3)
            # Points
            ax.scatter(data=xy_vals, x=resp_name, y=cyto, s='Cluster size', c='grey', edgecolor='black', zorder=3)
            #      Prediction intervall
            # ax.fill_between(x.transpose()[0], df_woweights['obs_ci_lower'], df_woweights['obs_ci_upper'], alpha=.5,
            #                 label='Prediction interval')

            # # Fit: The unweighted model
            ax.plot(x, ypred_woweights, ls='--', color='black', linewidth=2, label='Unweighted model')
            # Confidence intervall for 5 sigma: plotting error band
            ax.fill_between(x.transpose()[0], df_woweights['mean_ci_lower'].values,
                            df_woweights['mean_ci_upper'].values,
                            alpha=.1, label='5 - sigma interval', color='black', lw=0.1)
            #      Prediction intervall
            # ax.fill_between(x.transpose()[0], df_woweights['obs_ci_lower'], df_woweights['obs_ci_upper'], alpha=.5,
            #                 label='Prediction interval')

            # Axis params
            ax.set_xlabel(" ".join(["Responder Counts"]), fontsize=xy_fontsize)
            if cyto == 'IFNG':
                ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=xy_fontsize)
            else:
                ax.set_ylabel(" ".join([cyto, 'Counts']), fontsize=xy_fontsize)

            # set y-ticks
            if temp_df[cyto].max() < 10:
                ax.set_yticks(np.arange(0, temp_df[cyto].max() + 1, 1))
            elif (temp_df[cyto].max() >= 10) & (temp_df[cyto].max() < 50):
                ax.set_yticks(np.arange(0, temp_df[cyto].max() + 1, 5))
            elif (temp_df[cyto].max() >= 50) & (temp_df[cyto].max() < 100):
                ax.set_yticks(np.arange(0, temp_df[cyto].max() + 1, 10))
            elif (temp_df[cyto].max() >= 100) & (temp_df[cyto].max() < 200):
                ax.set_yticks(np.arange(0, temp_df[cyto].max() + 1, 40))
            elif (temp_df[cyto].max() >= 200) & (temp_df[cyto].max() < 600):
                ax.set_yticks(np.arange(0, temp_df[cyto].max() + 1, 50))
            else:
                ax.set_yticks(np.arange(0, temp_df[cyto].max() + 1, 100))

            # Add text: Correlation value and p-value
            ax.text(
                temp_df[resp_name].max() / 2 - temp_df[resp_name].max() / 10, temp_df[cyto].max(),
                'r = {:.2f}; p = {:.2e}'.format(corr_pval[0][2], corr_pval[0][3]),
                fontstyle='italic', fontsize=text_fontsize)
            ax.set_xlim([-0.5, temp_df[resp_name].max() + temp_df[resp_name].max() / 20])
            ax.set_ylim([-0.5, temp_df[cyto].max() + temp_df[cyto].max() / 20])

            # # Add text: slope and intercept term
            # ax.text(temp_df[resp_name].max() / 2 + temp_df[resp_name].max() / 5, 0,
            #         'm = {:.2e}; b = {:.2f}'.format(est_woweights_fitted.params[0], est_woweights_fitted.params[1]),
            #         fontstyle='italic', fontsize=text_fontsize)

            plt.tight_layout()
            # remove upper and right edge lines in plot
            sns.despine(ax=ax)

            # 3. Save figure
            fig.savefig(os.path.join(
                save_folder, "_".join(['Fig4A-C__weighted_transcripts_bulkapproach', cyto, resp_name, fileformat])))
            plt.close()

    # Save point size legend
    cluster_sizes = np.unique(pd.melt(df_counts_cytoresps, value_vars=list(genes_dict.keys()))['value'])
    plot_colorbar_legend.plot_standalone_legend(cluster_sizes, save_folder=save_folder, size_multiplier=size_multiplier)

    return dict_corr, df_counts_cytoresps


def main(save_folder, adata, corr_method):
    # parameter
    tissue_types = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']
    epidermis_layers = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']

    # 1. Get cytokines and responders
    t_cell_cytocines, cyto_resps_list, cytokine_responders = gene_lists.get_publication_cyto_resps()

    # 2. Subset adata to tissue_types of interest: upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS'
    if tissue_types:
        bool_col = adata.obs[tissue_types] == 1
        merged = bool_col.sum(axis=1)
        adata = adata[merged == 1]
        # Rename tissue region 'JUNCTION' to basal EPIDERMIS because some spots got both labels
        adata = ht.interface_to_epidermis(adata, tissue_layers='tissue_layer', epidermis_layers=epidermis_layers)

    # 3. Get counts of cyotkines and their responders in the EPIDERMIS
    df_correlation = get_bulk_cytoresp_counts(adata=adata, genes_dict=cytokine_responders, tissue_types=tissue_types)
    # Weighted by transcripts - goes into publication
    fig5__dict_weighted_transcripts_corr, df_counts_cytoresps = weighted_transcripts_correlation_plot(
        df_counts_cytoresps=df_correlation, genes_dict=cytokine_responders, save_folder=save_folder,
        corr_method=corr_method)

    return fig5__dict_weighted_transcripts_corr, df_counts_cytoresps


if __name__ == '__main__':
    today = date.today()

    # create saving folder in current project path
    savepath = os.path.join("..", "..", "..", "output", "Figure_4A", str(today))
    os.makedirs(savepath, exist_ok=True)

    # Load Raw anndata --> used for publication figure
    unpp_st_adata = sc.read(
        os.path.join("..", "..", "..", "adata_storage", "2022-04-08",
                     "Spatial Transcriptomics_unpp_cleaned_PsoADLP.h5"))

    main(save_folder=savepath, adata=unpp_st_adata, corr_method='spearman')
