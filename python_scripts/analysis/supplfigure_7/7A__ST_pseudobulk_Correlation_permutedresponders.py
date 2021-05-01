#!/usr/bin/env python
"""Calculate Pearson Correlation between Cytokines and their corresponding responder genes UMI-counts
    File name: 4A__ST_pseudobulk_aggregation_Correlation.py
    Author: Christina Hillig
    Date created: 23/11/2020
    Date last modified: 4/29/2021
    Python Version: 3.7
"""

from python_scripts.utils import gene_lists, add_observables

import scanpy as sc
import numpy as np
import pandas as pd
import os
from datetime import date

import pingouin

import seaborn as sns
import matplotlib.pyplot as plt


sc.set_figure_params(color_map='viridis')
fig_size = (8, 8)
xy_fontsize = 16
xy_ticks = 12
title_fontsize = 18
legend_fontsize = 14
text_fontsize = 14
fileformat = '.pdf'


def get_bulk_cytoresp_counts(adata, genes_dict, tissue_types=None):
    """Get per sample the counts of a cytokine and its responder genes

    Parameters
    ----------
    adata : annData
    genes_dict : dict
        keys are the cytokines and entries are the responder genes
    tissue_types : str

    Returns
    -------

    """

    if tissue_types:
        bool_col = adata.obs[tissue_types] == 1
        merged = bool_col.sum(axis=1)
        adata = adata[merged == 1]

    samples = np.unique(adata.obs['sample'])
    # get frequency
    col_names = list(genes_dict.keys())
    for gen in genes_dict.keys():
        col_names.append("_".join([gen, 'responder']))
    col_names.append('disease')

    df_corr = pd.DataFrame(index=samples, columns=col_names)

    for ind, sample in enumerate(samples):
        ad = adata[adata.obs['sample'] == sample]
        df_corr['disease'][ind] = np.unique(ad.obs['disease'].values)[0]
        for cytokine_gene in genes_dict:
            # 1. Get counts
            if "counts" in adata.layers.keys():
                counts_gene = ad.layers['counts'][:, np.where(ad.var.index == cytokine_gene)[0]]
            else:
                counts_gene = ad.X[:, np.where(ad.var.index == cytokine_gene)[0]]

            # sum counts of cytokine found in sample xy and add to dataframe
            counts_cyto_summed = np.asarray(counts_gene[:, 0], dtype=int).sum(axis=0)

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
            df_corr.at[sample, cytokine_gene] = counts_cyto_summed
            df_corr.at[sample, "_".join([cytokine_gene, 'responder'])] = counts_responders

    return df_corr


def correlation_plot(df_counts_cytoresps, genes_dict, save_folder):
    """
    Plot correlation of counts distribution for each cytokine over all samples

    :param df_counts_cytoresps:
    :param genes_dict:
    :param save_folder:
    :return:
    """

    # df_counts_cytoresps = df_counts_cytoresps.replace({'0': np.nan, 0: np.nan})
    df_counts_cytoresps = df_counts_cytoresps.replace({np.nan: 0})

    for cyto in genes_dict.keys():
        for cyto_reps in genes_dict.keys():
            resp_name = "_".join([cyto_reps, 'responder'])
            temp_df = df_counts_cytoresps[[cyto, resp_name, 'disease', 'marker']]

            # stats: hypothesis here will be that the counts of a cytokine is correlated to the counts of its responders
            sig_r = pingouin.corr(x=temp_df[resp_name], y=temp_df[cyto], method='pearson')

            # Plot Correlation
            #  2. Plot correlation
            fig, ax = plt.subplots(figsize=fig_size)
            ax.grid(False)
            sns.regplot(data=temp_df, x=resp_name, y=cyto, ax=ax, scatter=False, color="black", label=None)
            ax.scatter(data=temp_df, x=resp_name, y=cyto, c='k')

            # Axis params
            ax.set_xlabel(" ".join(["Responder Counts"]), fontsize=xy_fontsize)
            if cyto == 'IFNG':
                ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=xy_fontsize)
            else:
                ax.set_ylabel(" ".join([cyto, 'Counts']), fontsize=xy_fontsize)

            if temp_df.max()[0] < 10:
                ax.set_yticks(np.arange(0, temp_df.max()[0] + 1, 1))
                # Add text: Correlation value and p-value
                ax.text(temp_df.max()[1] / 2 - temp_df.max()[1] / 10, temp_df.max()[0],
                        'r = {:.2f}; p = {:.2e}'.format(sig_r['r'].values[0], sig_r['p-val'].values[0]),
                        fontstyle='italic', fontsize=text_fontsize)
            else:
                ax.set_yticks(np.arange(0, temp_df.max()[0] + 2, 2))
                # Add text: Correlation value and p-value
                ax.text(temp_df.max()[1] / 2 - temp_df.max()[1] / 10, temp_df.max()[0] + 1,
                        'r = {:.2f}; p = {:.2e}'.format(sig_r['r'].values[0], sig_r['p-val'].values[0]),
                        fontstyle='italic', fontsize=text_fontsize)
            ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 20])
            ax.set_ylim([-0.5, temp_df.max()[0] + temp_df.max()[0] / 20])

            plt.tight_layout()
            # remove upper and right edge lines in plot
            sns.despine(ax=ax)

            # 3. Save figure
            fig.savefig(os.path.join(save_folder, "_".join(['Fig4A', cyto, resp_name, fileformat])))
            plt.close()


def main(save_folder, adata, tissue_types):
    adata.obs['n_counts'] = adata.X.sum(1)
    # number of counts per spot in log
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    # number of genes per spot
    adata.obs['n_genes'] = (adata.X > 0).sum(1)
    mt_gene_mask = [gene.startswith('MT-') for gene in adata.var_names]
    adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1) / adata.obs['n_counts']

    # 1. Add meta data like which samples belong to which donor (optional)
    if "patient" not in adata.obs_keys():
        adata, tissue_cell_labels, disease_labels, lesion_labels = add_observables.add_metadata(adata)
        # 1.2 Remove spots having no tissue/cell labels (since 06.10.2020)
        adata = adata[np.where(adata.obs[tissue_cell_labels].to_numpy().any(axis=1))[0]]

    # 1. Get cytokines and responders
    _, cytokine_responders = gene_lists.get_permuted_respondergenes()

    # 2. Subset adata to tissue_types of interest: upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS'
    if tissue_types:
        bool_col = adata.obs[tissue_types] == 1
        merged = bool_col.sum(axis=1)
        adata = adata[merged == 1]

    # 3. Get counts of cyotkines and their responders in the EPIDERMIS
    df_correlation = get_bulk_cytoresp_counts(adata=adata, genes_dict=cytokine_responders, tissue_types=tissue_types)
    correlation_plot(df_counts_cytoresps=df_correlation, genes_dict=cytokine_responders, save_folder=save_folder)


if __name__ == '__main__':
    today = date.today()

    # parameter
    tissue_layers = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']

    # create saving folder in current project path
    savepath = os.path.join(os.environ['PYTHONPATH'].split(os.pathsep)[0], "output", "Figure_4A", str(today))
    os.makedirs(savepath, exist_ok=True)

    # Load Raw anndata --> used for publication figure
    unpp_st_adata = sc.read(os.path.join(savepath, "adata_storage/2020-10-06/st_adata_P15509_P16357_wo_4_7_unpp.h5"))

    main(save_folder=savepath, adata=unpp_st_adata, tissue_types=tissue_layers)
