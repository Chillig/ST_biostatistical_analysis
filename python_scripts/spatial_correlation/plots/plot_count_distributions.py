"""Visualise count distributions of responder and cytokine genes
    File name: plot_count_distributions.py
    Author: Christina Hillig
    Date created: May/01/2021
    Date last modified: May/01/2021
    Python Version: 3.7
"""

# import scripts
from python_scripts.spatial_correlation import corr_statistics as corr_stats
from python_scripts.spatial_correlation import helper_functions

# Plotting packages
import matplotlib.pyplot as plt
import seaborn as sns

# System specific
import os

# Calculation packages
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.stats as scstats


# Figure params
sc.set_figure_params(color_map='viridis')
fig_size, title_fontsize, axis_label_fontsize, legend_fontsize, fileformat, img_key, xy_ticks, text_fontsize = \
    helper_functions.figure_params()


def plot_counts(dfx, dfy, index_counter, index_counter_sample, cyto, sample, save_folder, distance=1):
    """Plot counts of cytokines and responders

    Parameters
    ----------
    dfx : pandas.Dataframe
    dfy : pandas.Dataframe
    index_counter : int
    index_counter_sample : int
    cyto : str
    sample : str
    save_folder : str
    distance : int

    Returns
    -------

    """
    fig, ax = plt.subplots(figsize=fig_size)
    ax.scatter(dfx[index_counter - index_counter_sample:index_counter],
               dfy[index_counter - index_counter_sample:index_counter], s=10, c='blue')
    ax.set_xlabel('responder genes', fontsize=axis_label_fontsize)
    ax.set_ylabel(" ".join([cyto, 'counts']), fontsize=axis_label_fontsize)
    ax.tick_params(labelsize=xy_ticks)
    ax.set_title(' '.join(['Distance', str(distance), sample, cyto]), fontsize=title_fontsize)
    plt.tight_layout()
    fig.savefig(
        os.path.join(save_folder, '_'.join(['Distance', str(distance), sample, cyto, "NotNormed", fileformat])))
    plt.close()


def plot_excluded_responder_spots(df_spot_counts, cytokines, save_folder):
    """Plot count distribution of excluded responder positives spots in the density clustering

    Parameters
    ----------
    df_spot_counts : pandas.Dataframe
    cytokines : list
    save_folder : str

    Returns
    -------

    """
    for ind, cyto in enumerate(cytokines):
        fig, ax = plt.subplots(figsize=fig_size)
        # Hide grid lines
        ax.grid(False)
        g = sns.distplot(df_spot_counts["_".join([cyto, 'responder'])], kde=False, color="xkcd:black", ax=ax)
        ax.set_xlabel(" ".join([cyto, 'responder counts']), fontsize=axis_label_fontsize)
        ax.set_ylabel('Count', fontsize=axis_label_fontsize)
        ax.set_title(" ".join([cyto, 'excluded responder counts']))
        g.set_yscale("log")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()
        fig.savefig(os.path.join(save_folder, "_".join(['Excluded_Responder_Distribution', cyto, fileformat])))
        plt.close()


def plot_included_responder_spots(df_spot_counts, cytokines, save_folder):
    """

    Parameters
    ----------
    df_spot_counts : pandas.Dataframe
    cytokines : list
    save_folder : str

    Returns
    -------

    """
    for ind, cyto in enumerate(cytokines):
        fig, ax = plt.subplots(figsize=fig_size)
        # Hide grid lines
        ax.grid(False)
        g = sns.distplot(df_spot_counts["_".join([cyto, 'responder'])], kde=False, color="xkcd:black", ax=ax)
        ax.set_xlabel(" ".join([cyto, 'responder counts']), fontsize=axis_label_fontsize)
        ax.set_ylabel('Count', fontsize=axis_label_fontsize)
        ax.set_title(" ".join([cyto, 'included responder counts']))
        g.set_yscale("log")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()
        fig.savefig(os.path.join(save_folder, "_".join(['Included_Responder_Distribution', cyto, fileformat])))
        plt.close()


def plot_distribution_respondercounts(adata, t_cell_cytocines, save_folder):
    """Plot count distribution of responder genes

    Parameters
    ----------
    adata : annData
    t_cell_cytocines : list
    save_folder : str

    Returns
    -------

    """
    for ind, cyto in enumerate(t_cell_cytocines):
        m_resp = adata.obs["_".join([cyto, 'Responders_counts'])] > 0
        fig, ax = plt.subplots(figsize=fig_size)
        # Hide grid lines
        ax.grid(False)
        g = sns.distplot(adata.obs["_".join([cyto, 'Responders_counts'])][m_resp], kde=False, color="xkcd:black", ax=ax)
        ax.set_xlabel(" ".join([cyto, 'responder counts']), fontsize=axis_label_fontsize)
        ax.set_ylabel('Count', fontsize=axis_label_fontsize)
        # ax.set_title("_".join([cyto, 'responder']))
        g.set_yscale("log")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()
        fig.savefig(os.path.join(save_folder, "_".join(['Responder_Distribution', cyto, fileformat])))
        plt.close()


def plot_responder_boxplot(df, cytokines, save_folder):
    """Plot boxplot of responder gene counts

    Parameters
    ----------
    df : pandas.Dataframe
    cytokines : list
    save_folder : str

    Returns
    -------

    """
    for cyto in cytokines:
        if cyto == 'IFNG':
            cyto_color = "#ff7f00"  # orange LICHEN
        elif cyto == "IL13":
            cyto_color = "#e41a1c"  # red AE
        else:
            cyto_color = "#377eb8"  # blue PSO
        col_resp_cytopos = "".join([cyto, '+ responder+ spot'])
        col_resp_cytoneg = "".join([cyto, '- responder+ spot'])

        # Apply test to check wether the responders are enriched in the cyto+ spots
        print(corr_stats.apply_wilcoxontest(df[col_resp_cytopos], df[col_resp_cytoneg]))

        # Two-sample Kolmogorov–Smirnov test
        ks_test = scstats.ks_2samp(df[col_resp_cytopos].values, df[col_resp_cytoneg].values)
        print(ks_test)

        df_melted = pd.melt(df[[col_resp_cytopos, col_resp_cytoneg]])

        fig, ax = plt.subplots(figsize=fig_size)
        # Turns off grid on the left Axis.
        ax.grid(False)
        sns.boxplot(x="variable", y="value", data=df_melted,
                    ax=ax, palette=[cyto_color, 'y']).set(xlabel='')

        ax.set_ylabel("Counts", fontsize=axis_label_fontsize)

        # statistical annotation
        x1, x2 = 0, 1
        yshift = df_melted['value'].max() / 15
        y, h, col = df_melted['value'].max() + yshift, yshift, 'k'
        plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        if ks_test.pvalue <= 0.05:
            plt.text((x1 + x2) * .5, y + h, "p = {:.2e}".format(ks_test.pvalue), ha='center', va='bottom', color=col)
        else:
            plt.text((x1 + x2) * .5, y + h, "ns", ha='center', va='bottom', color=col)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, "_".join([cyto, 'Boxplots_respondercounts.pdf'])))
        plt.close()


def plot_compare_inex_respondercounts(adata, df_includedcounts, df_excludedcounts, cytokines, save_folder):
    """Plot distribution of included (in) and excluded (ex) responder gene counts

    Parameters
    ----------
    adata : annData
    df_includedcounts : pandas.Dataframe
    df_excludedcounts : pandas.Dataframe
    cytokines : list
    save_folder : str

    Returns
    -------

    """
    for cyto in cytokines:
        if cyto == 'IFNG':
            cyto_color = "#ff7f00"  # orange LICHEN
        elif cyto == "IL13":
            cyto_color = "#e41a1c"  # red AE
        else:
            cyto_color = "#377eb8"  # blue PSO
        col_names = "_".join([cyto, 'responder'])

        # Add non-responder and non-cytokine count distribution
        # m_cytoresp = (adata.obs["_".join([cyto, 'Responders_counts'])].values > 0) & \
        #              (adata.obs["_".join([cyto, 'counts'])].values > 0)
        # counts = adata.obs['n_counts'][~m_cytoresp]

        # Apply test to check wether the responders are enriched in the cyto+ spots
        #  - only works if measurement have same samples size
        # print(corr_stats.apply_wilcoxontest(df_includedcounts[col_names].values, df_excludedcounts[col_names].values))

        # Two-sample Kolmogorov–Smirnov test
        ks_test = scstats.ks_2samp(df_includedcounts[col_names].values, df_excludedcounts[col_names].values)
        print(ks_test)

        included_colname = "included"
        excluded_colname = "excluded"

        df_includedcounts = df_includedcounts.rename(columns={col_names: included_colname})
        df_excludedcounts = df_excludedcounts.rename(columns={col_names: excluded_colname})
        # df_otherscounts = pd.DataFrame(counts.values, columns=['Others'])

        # df_merged = pd.concat([df_includedcounts[included_colname],
        #                        df_excludedcounts[excluded_colname], df_otherscounts['Others']], axis=1)
        df_merged = pd.concat([df_includedcounts[included_colname],
                               df_excludedcounts[excluded_colname]], axis=1)
        df_melted = pd.melt(df_merged)

        fig, ax = plt.subplots(figsize=fig_size)
        # Turns off grid on the left Axis.
        ax.grid(False)
        sns.boxplot(x="variable", y="value", data=df_melted, ax=ax, palette=[cyto_color, 'y', 'gray']).set(xlabel='')
        ax.set_ylabel("Counts", fontsize=axis_label_fontsize)

        # statistical annotation
        x1, x2 = 0, 1
        yshift = df_melted['value'].max() / 15
        y, h, col = df_melted['value'].max() + yshift, yshift, 'k'
        plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        if ks_test.pvalue <= 0.05:
            plt.text((x1 + x2) * .5, y + h, "p = {:.2e}".format(ks_test.pvalue), ha='center', va='bottom', color=col)
        else:
            plt.text((x1 + x2) * .5, y + h, "ns", ha='center', va='bottom', color=col)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        plt.title('{} responder counts'.format(cyto))

        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, "_".join([cyto, 'Boxplots_InExcluded_respondercounts.pdf'])))
        plt.close()


def plot_responder_distributions(adata, df_included_responder, df_excluded_responder, t_cell_cytocines, save_folder):
    """Distribution plot of included (in conditional density cluster) and excluded Responder counts

    Parameters
    ----------
    adata : annData
    df_included_responder : pandas.Dataframe
    df_excluded_responder : pandas.Dataframe
    t_cell_cytocines : list
    save_folder : str

    Returns
    -------

    """
    for ind, cyto in enumerate(t_cell_cytocines):
        col_names = "_".join([cyto, 'responder'])
        m_resp = adata.obs["_".join([cyto, 'Responders_counts'])] > 0
        fig, ax = plt.subplots(figsize=fig_size)
        # Hide grid lines
        ax.grid(False)
        bins = np.linspace(adata.obs["_".join([cyto, 'Responders_counts'])][m_resp].min(),
                           adata.obs["_".join([cyto, 'Responders_counts'])][m_resp].max(), 40)
        sns.distplot(adata.obs["_".join([cyto, 'Responders_counts'])][m_resp], kde=False,
                     color="k", ax=ax, label="All responder counts", hist_kws={"alpha": 0.7}, bins=bins)
        sns.distplot(df_excluded_responder[col_names], kde=False, color="darkorange",
                     ax=ax, label="Excluded responder counts", hist_kws={"alpha": 0.7}, bins=bins)
        g = sns.distplot(df_included_responder[col_names], kde=False, color="red", ax=ax,
                         label="Included responder counts", hist_kws={"alpha": 0.7}, bins=bins)
        ax.set_xlabel(" ".join([cyto, 'responder counts']), fontsize=axis_label_fontsize)
        ax.set_ylabel('Count', fontsize=axis_label_fontsize)
        g.set_yscale("log")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.legend()

        plt.tight_layout()
        fig.savefig(os.path.join(save_folder, "_".join(['Compare_Responder_Distribution', cyto, fileformat])))
        plt.close()
