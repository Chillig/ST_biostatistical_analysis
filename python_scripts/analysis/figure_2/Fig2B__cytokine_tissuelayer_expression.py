#!/usr/bin/env python
"""Get Distribution of UMI-counts of cytokines
    File name: Fig2B__cytokine_tissuelayer_expression.py
    Author: Christina Hillig
    Date created: 23/11/2020
    Date last modified: 3/21/2021
    Python Version: 3.7

    Script description:
    1. Get counts of cytokine in tissue sub-groups of EPIDERMIS and DERMIS
    2. Apply Wilcoxon signed-rank test to test if cytokines are enrichment in specific tissue parts
"""
from python_scripts.utils import gene_lists

import scanpy as sc
import os
from datetime import date
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from scipy import stats


# figure properties
figure_size = (8, 8)
xy_fontsize = 14
xy_ticks = 12
title_fontsize = 16
legend_fontsize = 10
fileformat = '.pdf'


def get_dotsize(freq):
    """Get the size of the dots dependent on the frequency value

    Parameters
    ----------
    freq : numpy.array

    Returns
    -------

    """
    dotsize_range = np.zeros_like(freq)
    dotsize_range[(freq >= 1) & (freq < 10)] = 10
    dotsize_range[(freq >= 10) & (freq < 100)] = 100
    dotsize_range[(freq >= 100) & (freq < 1000)] = 400
    dotsize_range[(freq >= 1000)] = 800

    return dotsize_range, [0, 1, 10, 100, 1000]


def plot_umicount_dotplot(df, signatures_hkg, title, ylabels, xlabels, save_folder, filename):
    """Plot UMI-counts distribution of cytokine and housekeeping gene, GAPDH, in a Dotplot

    Parameters
    ----------
    df
    signatures_hkg
    title
    ylabels
    xlabels
    save_folder
    filename

    Returns
    -------

    """
    fig, ax = plt.subplots(nrows=1, ncols=1, facecolor='w', edgecolor='k', figsize=figure_size)
    fig.subplots_adjust(hspace=0.8, wspace=0.25)
    ax.grid(False)
    ax.set_title(title, fontsize=title_fontsize)
    ax.scatter('x-axis', 'y-axis', data=df, s='dot_size', c='color', linewidths=.5, edgecolors='k')
    ax.set_xlim([-0.2, len(signatures_hkg) - 0.8])
    ax.set_xticks(np.arange(0, len(signatures_hkg), 1.0))
    ax.set_yticks(np.arange(0, len(ylabels), 1.0))
    ax.set_yticklabels(ylabels, fontsize=xy_fontsize)
    ax.set_xticklabels(xlabels, fontsize=xy_fontsize)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend:
    # produce a legend with a cross section of sizes from the scatter
    # make a legend:
    fake_dotsize, rpw = get_dotsize(np.arange(1e4))
    for ind, pw in enumerate(np.unique(fake_dotsize)):
        if pw != 0:
            plt.scatter([], [], s=pw, c="k", label=" >{}".format(rpw[ind]))

    handle, leg = plt.gca().get_legend_handles_labels()
    plt.legend(handle[1:], leg[1:], labelspacing=2.5, title="UMI-counts",
               title_fontsize=title_fontsize,
               fontsize=legend_fontsize,  # size of title
               bbox_to_anchor=(1.001, 0.8,), bbox_transform=fig.transFigure,
               borderpad=1, scatterpoints=1, handletextpad=1, handlelength=1,
               frameon=True, framealpha=0.6, edgecolor="k", facecolor="w")

    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "".join([filename, fileformat])), bbox_inches='tight')
    plt.close()


def plot_umicounts_count(df, title, sample_name, save_folder, ylabels, hstacks):
    """Plot count of UMI-counts

    Parameters
    ----------
    df : pandas.Dataframe
    title : str
    sample_name : str
    save_folder : str
    ylabels : list
    hstacks : list

    Returns
    -------

    """
    # plot count of cytokines
    fig, ax = plt.subplots(figsize=figure_size)
    # We need to draw the canvas, otherwise the labels won't be positioned and
    # won't have values yet.
    fig.canvas.draw()
    height = 0.2
    operation = [np.subtract, np.add, np.add, np.subtract]
    shift = [1.55, 1.45, 0.45, 0.55]
    np.mod(2, 3)
    for ind, gene in enumerate(hstacks):
        m_gene = df['Genes'] == gene
        yvalue = operation[ind](df['y-axis'][m_gene], np.multiply(shift[ind], height))
        ax.barh(operation[ind](df['y-axis'][m_gene], np.multiply(shift[ind], height)), align='center',
                height=height, width=df['frequency'][m_gene], color=df['color'][m_gene], label=gene)
        # Annotate stacks
        for ind_y, v in enumerate(df['y-axis'][m_gene]):
            if df['frequency'][m_gene].iloc[ind_y] != 0:
                ax.text(df['frequency'][m_gene].iloc[ind_y] + 0.1, yvalue.iloc[ind_y] + 0.05,
                        str(df['frequency'][m_gene].iloc[ind_y]), color='k', fontweight='bold', fontsize=6)
    ax.set_xscale('log')
    ax.set_xlabel("Counts", labelpad=14)
    ax.set_ylabel("Tissue layers", labelpad=14)
    ax.set_yticks(np.arange(0, len(list(np.unique(df['tissue_types'])))))
    ax.set_yticklabels(ylabels)
    ax.invert_yaxis()
    # remove borders from plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    # set title
    plt.title(" ".join(["Counts of", title]), y=1.02)
    plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "_".join([sample_name, "Count", fileformat])))
    plt.close()


def plot_signature_umicounts(df_l, df_nl, save_folder):
    """Plot UMI-counts distribution of signature cytokines

    Parameters
    ----------
    df_l : pandas.Dataframe
    df_nl : pandas.Dataframe
    save_folder : str

    Returns
    -------

    """

    # 1. Select tissue types of interest
    tissue_types = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS',
                    'DERdepth1', 'DERdepth2', 'DERdepth3', 'DERdepth4', 'DERdepth5', 'DERdepth6', 'DERdepth7']

    signatures = gene_lists.get_color_signaturegenes()
    signatures_list_l = list(signatures.keys())
    df_freq_l = pd.DataFrame(columns=['x-axis', 'y-axis', "Genes", 'tissue_types', 'color', 'frequency',
                                      'log_frequency'])
    for ind, tissue in enumerate(tissue_types):
        freq = df_l.loc[tissue]
        dotsize_range, _ = get_dotsize(freq=freq)
        df_long = pd.DataFrame({'x-axis': np.arange(0, len(signatures_list_l), 1).astype(np.int64),
                                'y-axis': np.ones(len(signatures_list_l)).astype(np.int64) * ind,
                                "Genes": signatures_list_l,
                                'tissue_types': [tissue] * len(signatures_list_l),
                                'color': list(signatures.values()),
                                'frequency': freq.astype(np.int64), 'log_frequency': np.log10(freq.astype(np.float64)),
                                'dot_size': dotsize_range.astype(np.int64)})
        df_freq_l = df_freq_l.append(df_long, ignore_index=True)

    signatures_list_nl = list(signatures.keys())
    df_freq_nl = pd.DataFrame(columns=['x-axis', 'y-axis', "Genes", 'tissue_types', 'color', 'frequency',
                                       'log_frequency'])
    for ind, tissue in enumerate(tissue_types):
        freq = df_nl.loc[tissue]
        dotsize_range, _ = get_dotsize(freq=freq.values)
        df_long = pd.DataFrame({'x-axis': np.arange(0, len(signatures_list_nl), 1).astype(np.int64),
                                'y-axis': np.ones(len(signatures_list_nl)).astype(np.int64) * ind,
                                "Genes": signatures_list_nl,
                                'tissue_types': [tissue] * len(signatures_list_nl),
                                'color': list(signatures.values()),
                                'frequency': freq.astype(np.int64), 'log_frequency': np.log10(freq.astype(np.float64)),
                                'dot_size': dotsize_range.astype(np.int64)})
        df_freq_nl = df_freq_nl.append(df_long, ignore_index=True)

    """ Figure 2B: Dotplot  """
    # Lesioned Skin test plot
    plot_umicount_dotplot(df=df_freq_l, signatures_hkg=signatures_list_l,
                          title="Lesioned Skin", ylabels=tissue_types, xlabels=list(signatures.keys()),
                          save_folder=save_folder, filename='Figure2a__L_Skin_dotplot')

    # Non Lesioned Skin test plot
    plot_umicount_dotplot(df=df_freq_nl, signatures_hkg=signatures_list_nl,
                          title="Non Lesioned Skin", ylabels=tissue_types, xlabels=list(signatures.keys()),
                          save_folder=save_folder, filename='Figure2a__NL_Skin_dotplot')


def divide_lnl_adata(adata, save_folder):
    """Read out counts per tissue and separate them into lesion and non lesion skin

    Parameters
    ----------
    adata : annData
    save_folder : str

    Returns
    -------

    """
    # 1.1 signature cytokines for diseases LICHEN, AE and PSO and the Housekeeping gene GAPDH
    sig_types = ['IFNG', 'IL13', 'IL17A', 'GAPDH']

    # 1.2 Select tissue types of interest
    tissue_types = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS',
                    'DERdepth1', 'DERdepth2', 'DERdepth3', 'DERdepth4', 'DERdepth5', 'DERdepth6', 'DERdepth7']
    # Get available tissue types
    available_tissues = list(set(adata.obs_keys()) & set(tissue_types))
    available_tissues = sorted(available_tissues)

    # divide into lesion and non lesional samples
    m_nonlesion = (adata.obs['NON LESIONAL'] == 1).values
    m_lesion = (adata.obs['NON LESIONAL'] == 0).values

    # create dataframes
    df_nl = pd.DataFrame(index=available_tissues, columns=list(sig_types))
    df_l = pd.DataFrame(index=available_tissues, columns=list(sig_types))

    # 2.1 Get boolean of tissue parts ==1 in adata.obs => bool matrix of barcode/spot if belongs to tissue part
    m_barcodes_tissues = (adata.obs[available_tissues] == 1).values

    info_lesionedskin = pd.DataFrame(index=available_tissues, columns=list(sig_types))
    info_non_lesionedskin = pd.DataFrame(index=available_tissues, columns=list(sig_types))

    for gene in sig_types:
        # 3. Get signature spot positions
        if gene in adata.var_names:
            if "counts" in adata.layers.keys():
                gene_matrix = np.copy(adata.layers['counts'])[:, np.where(adata.var.index == gene)[0]]
            else:
                gene_matrix = np.copy(adata.X)[:, np.where(adata.var.index == gene)[0]]

            for ind, tissue_part in enumerate(available_tissues):
                obs_name = "_".join([gene, 'counts', tissue_part])
                # get counts of gene in each tissue type
                counts_nonlesion_tissue = gene_matrix[:, 0][m_barcodes_tissues[:, ind] * m_nonlesion]
                counts_lesion_tissue = gene_matrix[:, 0][m_barcodes_tissues[:, ind] * m_lesion]
                counts_tissue = gene_matrix[:, 0][m_barcodes_tissues[:, ind]]

                # save to adata and data frame
                df_l[gene].loc[tissue_part] = counts_lesion_tissue.sum()
                df_nl[gene].loc[tissue_part] = counts_nonlesion_tissue.sum()

                info_lesionedskin[gene].loc[tissue_part] = np.count_nonzero(counts_lesion_tissue)
                info_non_lesionedskin[gene].loc[tissue_part] = np.count_nonzero(counts_nonlesion_tissue)

                adata.obs[obs_name] = 0
                adata.obs[obs_name].loc[m_barcodes_tissues[:, ind]] = counts_tissue

    # save data frame with information
    info_non_lesionedskin.to_csv(os.path.join(save_folder, "info_non-lesionedSkin.csv"))
    info_lesionedskin.to_csv(os.path.join(save_folder, "info_lesionedSkin.csv"))
    df_l.to_csv(os.path.join(save_folder, "counts_lesionedSkin.csv"))
    df_nl.to_csv(os.path.join(save_folder, "counts_non-lesionedSkin.csv"))

    return adata, df_l, df_nl


def apply_wilcoxontest(df_highcounts, df_lowcounts, save_folder):
    """Apply Wilcoxon signed-rank test to check if the two counts distributions are significantly different

    Parameters
    ----------
    df_highcounts : pandas.Dataframe
    df_lowcounts : pandas.Dataframe
    save_folder : str

    Returns
    -------

    """
    stats.probplot(df_highcounts, dist="norm", plot=plt)
    plt.title("High Counts Q-Q Plot")
    plt.savefig(os.path.join(save_folder, 'QQ-plot_High_Counts.png'))

    stats.probplot(df_lowcounts, dist="norm", plot=plt)
    plt.title("Lower Counts Q-Q Plot")
    plt.savefig(os.path.join(save_folder, 'QQ-plot_Lower_Counts.png'))

    # check if data is normal distributed
    w, ps_hc = stats.shapiro(df_highcounts)
    w, ps_lc = stats.shapiro(df_lowcounts)
    # if p-value < 0.05 -> variable violates the assumption of normality => Use Wilcoxon signed rank test
    if (ps_hc < 0.05) & (ps_lc < 0.05):
        t, p_w = stats.wilcoxon(df_highcounts, df_lowcounts)
        print(p_w)
        return p_w
    else:
        print("Distributions are not significantly different")


def main(adata, save_folder):
    """

    Parameters
    ----------
    adata : annData
    save_folder : str

    Returns
    -------

    """
    # 1. Split data set into lesional and non-lesional
    adata, df_l, df_nl = divide_lnl_adata(adata, save_folder=save_folder)

    """Paper Figure 2B: Test if cytokines are enriched in any tissue parts"""
    plot_signature_umicounts(df_l, df_nl, save_folder)

    # mask for lesional samples
    m_lesion = [adata.obs['NON LESIONAL'] == 0]

    # Read out cytokine counts in certain tissue layers and norm them by the number of involved tissue types
    df_hc_ifng = (adata.obs['IFNG_counts_basal EPIDERMIS'].values[tuple(m_lesion)] +
                  adata.obs['IFNG_counts_DERdepth1'].values[tuple(m_lesion)])
    df_lc_ifng = (adata.obs['IFNG_counts_upper EPIDERMIS'].values[tuple(m_lesion)] +
                  adata.obs['IFNG_counts_middle EPIDERMIS'].values[tuple(m_lesion)] +
                  adata.obs['IFNG_counts_DERdepth2'].values[tuple(m_lesion)] +
                  adata.obs['IFNG_counts_DERdepth3'].values[tuple(m_lesion)] +
                  adata.obs['IFNG_counts_DERdepth4'].values[tuple(m_lesion)] +
                  adata.obs['IFNG_counts_DERdepth5'].values[tuple(m_lesion)] +
                  adata.obs['IFNG_counts_DERdepth7'].values[tuple(m_lesion)])

    # Apply Wilcoxon test to check if IFN-g is enriched in specific tissue layers
    apply_wilcoxontest(df_highcounts=df_hc_ifng / 2, df_lowcounts=df_lc_ifng / 6, save_folder=save_folder)

    df_hc_il13 = (adata.obs['IL13_counts_middle EPIDERMIS'].values[tuple(m_lesion)] +
                  adata.obs['IL13_counts_basal EPIDERMIS'].values[tuple(m_lesion)] +
                  adata.obs['IL13_counts_DERdepth1'].values[tuple(m_lesion)] +
                  adata.obs['IL13_counts_DERdepth2'].values[tuple(m_lesion)]) / 4
    df_lc_il13 = (adata.obs['IL13_counts_DERdepth3'].values[tuple(m_lesion)] +
                  adata.obs['IL13_counts_DERdepth4'].values[tuple(m_lesion)] +
                  adata.obs['IL13_counts_DERdepth5'].values[tuple(m_lesion)] +
                  adata.obs['IL13_counts_DERdepth6'].values[tuple(m_lesion)] +
                  adata.obs['IL13_counts_DERdepth7'].values[tuple(m_lesion)] +
                  adata.obs['IL13_counts_upper EPIDERMIS'].values[tuple(m_lesion)]) / 6

    # Apply Wilcoxon test to check if IL-13 is enriched in specific tissue layers
    apply_wilcoxontest(df_highcounts=df_hc_il13, df_lowcounts=df_lc_il13, save_folder=save_folder)

    df_hc_il17a = (adata.obs['IL17A_counts_basal EPIDERMIS'].values[tuple(m_lesion)] +
                   adata.obs['IL17A_counts_upper EPIDERMIS'].values[tuple(m_lesion)] +
                   adata.obs['IL17A_counts_middle EPIDERMIS'].values[tuple(m_lesion)] +
                   adata.obs['IL17A_counts_DERdepth1'].values[tuple(m_lesion)]) / 4
    df_lc_il17a = (adata.obs['IL17A_counts_DERdepth3'].values[tuple(m_lesion)] +
                   adata.obs['IL17A_counts_DERdepth4'].values[tuple(m_lesion)] +
                   adata.obs['IL17A_counts_DERdepth5'].values[tuple(m_lesion)] +
                   adata.obs['IL17A_counts_DERdepth6'].values[tuple(m_lesion)] +
                   adata.obs['IL17A_counts_DERdepth7'].values[tuple(m_lesion)] +
                   adata.obs['IL17A_counts_DERdepth2'].values[tuple(m_lesion)]) / 7

    # Apply Wilcoxon test to check if IL-17A is enriched in specific tissue layers
    apply_wilcoxontest(df_highcounts=df_hc_il17a, df_lowcounts=df_lc_il17a, save_folder=save_folder)


if __name__ == '__main__':
    today = date.today()
    # create saving folder
    save_path = os.path.join("..", "..", "..", "output", "Figure_2B", str(today))
    os.makedirs(save_path, exist_ok=True)

    # Load unpreprocessed but not normalised annData object
    unpp_adata = sc.read(
        os.path.join("..", "..", "..", "adata_storage", "2022-04-08",
                     "Spatial Transcriptomics_unpp_cleaned_PsoADLP.h5"))

    main(adata=unpp_adata, save_folder=save_path)
