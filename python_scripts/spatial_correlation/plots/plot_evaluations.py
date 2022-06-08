"""Evaluation metrics
    File name: plot_evaluations.py
    Author: Christina Hillig
    Date created: May/01/2021
    Date last modified: May/01/2021
    Python Version: 3.7
"""

# import scripts
from python_scripts.spatial_correlation import helper_functions
from python_scripts.utils import gene_lists

# Plotting packages
import matplotlib.pyplot as plt
import seaborn as sns

# System specific
import os

# Calculation packages
import scanpy as sc
import numpy as np
import pandas as pd


# Figure params
sc.set_figure_params(color_map='viridis')
fig_size, title_fontsize, axis_label_fontsize, legend_fontsize, fileformat, img_key, xy_ticks, text_fontsize = \
    helper_functions.figure_params()


def plot_evaluate_distance(significance, cytokines, min_radius, save_folder):
    """Elbow plot for best distance/ radius evaluation

    Parameters
    ----------
    significance : list
    cytokines : list
    min_radius: int
    save_folder : str

    Returns
    -------

    """
    # Evaluate distance via elbow plot
    significance = np.array(significance).T
    sig_pearson = []
    sig_spearman = []
    for val in significance:
        sig_pearson.append(val['pearson'])
        sig_spearman.append(val['spearman'])

    sig_pearson = np.array(sig_pearson)

    sig_spearman = np.array(sig_spearman)

    # load cytokine to color
    cyto_color = helper_functions.get_color_signaturegenes()

    if min_radius > 0:
        x_vals = np.arange(min_radius, sig_pearson.shape[0] + 1)
    else:
        x_vals = np.arange(min_radius, sig_pearson.shape[0])

    if sig_pearson.T[1:].astype('float').min() < 0:
        ymin_p = -1
    else:
        ymin_p = 0

    if sig_spearman.T[1:].astype('float').min() < 0:
        ymin_s = -1
    else:
        ymin_s = 0

    fig_pval, ax_pval = plt.subplots(figsize=fig_size)
    # ax_pval.grid(False)
    fig_corr, ax_corr = plt.subplots(figsize=fig_size)
    ax_corr.set_ylim([ymin_p, 1])
    # ax_corr.grid(False)
    for ind, cyto in enumerate(cytokines):  # KeyError: 'IFNG_IL13_responder'
        mask = sig_pearson.T[:, ind, :].T[:, 2].astype('float') < 0.05
        ind_notsigpval = np.where(mask == False)[0]
        ind_sigpval = np.where(mask == True)[0]

        ax_pval.plot(x_vals, -np.log10(sig_pearson.T[:, ind, :].T[:, 2].astype('float')),
                     linestyle='-', c=cyto_color[cyto], label=cyto)
        # Highlight significant markers with triangle and non significant ones with unfilled circle
        if len(ind_sigpval) > 0:
            ax_pval.scatter(x_vals[mask], -np.log10(sig_pearson.T[:, ind, :].T[:, 2].astype('float'))[mask],
                            linestyle='-', marker='^', c=cyto_color[cyto])
        if len(ind_notsigpval) > 0:
            ax_pval.scatter(x_vals[~mask], -np.log10(sig_pearson.T[:, ind, :].T[:, 2].astype('float'))[~mask],
                            linestyle='-', marker='o', c=cyto_color[cyto], facecolors='none')
        ax_pval.set_xlabel('Radius', fontsize=axis_label_fontsize)
        ax_pval.set_ylabel(r'-log$_{10}$(p-values)', fontsize=axis_label_fontsize)
        ax_pval.set_xticks(x_vals)
        sns.despine(ax=ax_pval)

        # Plot correlation vs radius
        ax_corr.plot(x_vals, sig_pearson.T[:, ind, :].T[:, 1].astype('float'),
                     linestyle='-', c=cyto_color[cyto], label=cyto)
        # Highlight significant markers with triangle and non significant ones with unfilled circle
        if len(ind_sigpval) > 0:
            ax_corr.scatter(x_vals[mask], sig_pearson.T[:, ind, :].T[:, 1].astype('float')[mask],
                            linestyle='-', marker='^', c=cyto_color[cyto])
        if len(ind_notsigpval) > 0:
            ax_corr.scatter(x_vals[~mask], sig_pearson.T[:, ind, :].T[:, 1].astype('float')[~mask],
                            linestyle='-', marker='o', c=cyto_color[cyto], facecolor='white')
        ax_corr.set_xlabel('Radius', fontsize=axis_label_fontsize)
        ax_corr.set_ylabel(r'Correlation value', fontsize=axis_label_fontsize)
        ax_corr.set_xticks(x_vals)
        sns.despine(ax=ax_corr)

    # ax_pval.legend()
    # fig_pval.savefig(os.path.join(save_folder, "_".join(['PearsonPval_vs_Radius_Evaluation_wogrid', fileformat])))
    # plt.close(fig=fig_pval)
    # ax_corr.legend()
    # fig_corr.savefig(os.path.join(save_folder, "_".join(['PearsonCorr_vs_Radius_Evaluation_wogrid', fileformat])))
    # plt.close(fig=fig_corr)

    leg = ax_pval.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, fancybox=True, shadow=False)
    fig_pval.savefig(os.path.join(save_folder, "_".join(['PearsonPval_vs_Radius_Evaluation', fileformat])),
                     bbox_inches='tight',  bbox_extra_artists=(leg,))
    plt.close(fig=fig_pval)
    leg = ax_corr.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, fancybox=True, shadow=False)
    fig_corr.savefig(os.path.join(save_folder, "_".join(['PearsonCorr_vs_Radius_Evaluation', fileformat])),
                     bbox_inches='tight',  bbox_extra_artists=(leg,))
    plt.close(fig=fig_corr)

    fig_pval, ax_pval = plt.subplots(figsize=fig_size)
    # ax_pval.grid(False)
    fig_corr, ax_corr = plt.subplots(figsize=fig_size)
    ax_corr.set_ylim([ymin_s, 1])
    # ax_corr.grid(False)
    for ind, cyto in enumerate(cytokines):
        mask = sig_spearman.T[:, ind, :].T[:, 2].astype('float') < 0.05
        ind_notsigpval = np.where(mask == False)[0]
        ind_sigpval = np.where(mask == True)[0]

        ax_pval.plot(x_vals, -np.log10(sig_spearman.T[:, ind, :].T[:, 2].astype('float')),
                     linestyle='-', c=cyto_color[cyto], label=cyto)
        # Highlight significant markers with triangle and non significant ones with unfilled circle
        if len(ind_sigpval) > 0:
            ax_pval.scatter(x_vals[mask], -np.log10(sig_spearman.T[:, ind, :].T[:, 2].astype('float'))[mask],
                            linestyle='-', marker='^', c=cyto_color[cyto])
        if len(ind_notsigpval) > 0:
            ax_pval.scatter(x_vals[~mask], -np.log10(sig_spearman.T[:, ind, :].T[:, 2].astype('float'))[~mask],
                            linestyle='-', marker='o', c=cyto_color[cyto], facecolors='none')
        ax_pval.set_xlabel('Radius', fontsize=axis_label_fontsize)
        ax_pval.set_ylabel(r'-log$_{10}$(p-values)', fontsize=axis_label_fontsize)
        ax_pval.set_xticks(x_vals)
        sns.despine(ax=ax_pval)

        # Plot correlation vs radius
        ax_corr.plot(x_vals, sig_spearman.T[:, ind, :].T[:, 1].astype('float'),
                     linestyle='-', c=cyto_color[cyto], label=cyto)
        # Highlight significant markers with triangle and non significant ones with unfilled circle
        if len(ind_sigpval) > 0:
            ax_corr.scatter(x_vals[mask], sig_spearman.T[:, ind, :].T[:, 1].astype('float')[mask],
                            linestyle='-', marker='^', c=cyto_color[cyto])
        if len(ind_notsigpval) > 0:
            ax_corr.scatter(x_vals[~mask], sig_spearman.T[:, ind, :].T[:, 1].astype('float')[~mask],
                            linestyle='-', marker='o', c=cyto_color[cyto], facecolor='white')
        ax_corr.set_xlabel('Radius', fontsize=axis_label_fontsize)
        ax_corr.set_ylabel(r'Correlation value', fontsize=axis_label_fontsize)
        ax_corr.set_xticks(x_vals)
        sns.despine(ax=ax_corr)

    # ax_pval.legend()
    # fig_pval.savefig(os.path.join(save_folder, "_".join(['SpearmanPval_vs_Radius_Evaluation_wogrid', fileformat])))
    # plt.close(fig=fig_pval)
    # ax_corr.legend()
    # fig_corr.savefig(os.path.join(save_folder, "_".join(['SpearmanCorr_vs_Radius_Evaluation_wogrid', fileformat])))
    # plt.close(fig=fig_corr)

    leg = ax_pval.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, fancybox=True, shadow=False)
    fig_pval.savefig(os.path.join(save_folder, "_".join(['SpearmanPval_vs_Radius_Evaluation', fileformat])),
                     bbox_inches='tight',  bbox_extra_artists=(leg,))
    plt.close(fig=fig_pval)
    leg = ax_corr.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, fancybox=True, shadow=False)
    fig_corr.savefig(os.path.join(save_folder, "_".join(['SpearmanCorr_vs_Radius_Evaluation', fileformat])),
                     bbox_inches='tight',  bbox_extra_artists=(leg,))
    plt.close(fig=fig_corr)


def plot_responder_vs_radius(counts_dict: dict, conditionalgenes_responders: dict, radii: list, save_folder: str):
    # Normalise Responder counts by by number of clusters
    df_respcounts_radius = pd.DataFrame(columns=['normed_responder', 'radius', 'cytokine'])
    for radius in radii:
        for cyto in conditionalgenes_responders.keys():
            df_temp = pd.DataFrame(columns=['normed_responder', 'radius', 'cytokine'])
            for ind, specimen in enumerate(counts_dict[radius]['Specimen'].unique()):
                df_temp_specimen = pd.DataFrame(columns=['normed_responder', 'radius', 'cytokine', 'num_clusters'])
                df_specimen = counts_dict[radius][counts_dict[radius].loc[:, 'Specimen'] == specimen].copy()
                # Number of clusters for a specific cytokine on a specimen
                num_clusters = df_specimen[~df_specimen['{}_responder'.format(cyto)].isna()].shape[0]
                numspots_cluster = df_specimen[~df_specimen['{}_responder'.format(cyto)].isna()][
                    'Cluster_num_spots']

                # Responder counts normed by clusters
                if np.count_nonzero(numspots_cluster) > 0:
                    numspots_cluster = numspots_cluster.replace(0, np.nan)
                    normed_resp_counts = df_specimen['{}_responder'.format(cyto)].astype(float) / numspots_cluster
                    normed_resp_counts = normed_resp_counts.sum()
                else:
                    normed_resp_counts = np.nan

                df_temp_specimen.loc[ind, 'normed_responder'] = normed_resp_counts
                df_temp_specimen.loc[ind, 'radius'] = radius
                df_temp_specimen.loc[ind, 'cytokine'] = cyto
                df_temp_specimen.loc[ind, 'num_clusters'] = num_clusters

                df_temp = pd.concat([df_temp, df_temp_specimen])

            df_respcounts_radius = pd.concat([df_respcounts_radius, df_temp], axis=0)

    # Draw radius vs Responder counts normed by clusters
    for cyto in conditionalgenes_responders.keys():
        fig, ax = plt.subplots()
        ax.grid(False)
        sns.pointplot(x="radius", y="normed_responder", ci="sd", capsize=0.1,
                      data=df_respcounts_radius[df_respcounts_radius['cytokine'] == cyto],
                      dodge=True, join=False, ax=ax)
        sns.despine(ax=ax)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Normed {} responder counts'.format(cyto.split('_')[0]))
        plt.tight_layout()
        fig.savefig(os.path.join(save_folder, '{}__Radius_vs_normed_bynspots_Respcounts.pdf'.format(cyto)))
        plt.close(fig=fig)
