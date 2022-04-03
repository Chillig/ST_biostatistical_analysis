#!/usr/bin/env python
"""Signed log10-transformed p-value plot comparing DGE Analysis results of ST and scRNA-seq data sets
    File name: Fig3E__STSC_Correlation_IL17A.py
    Author: Christina Hillig
    Date created: March/06/2021
    Date last modified: April/29/2021
    Python Version: 3.7
"""

from python_scripts.utils import gene_lists

from datetime import date
import os
import pandas as pd
import numpy as np
import pingouin

import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.patches import Rectangle
import plotly.graph_objects as go

dotsize = 9
fig_size = (8, 8)
xy_fontsize = 16
xy_ticks = 12
title_fontsize = 18
legend_fontsize = 14
text_fontsize = 12
file_format = '.pdf'

# source: https://www.genomics-online.com/resources/16/5049/housekeeping-genes/
_house_keeping_genes = ['ACTB', 'GAPDH', 'PGK1', 'PPIA', 'RPLP0', 'ARBP', 'B2M', 'YWHAZ', 'SDHA', 'TFRC', 'GUSB',
                        'HMBS', 'HPRT1', 'TBP']


def check_masks(df_st, df_sc, value, sig_13_mask, sig_24_mask, not_sig_mask, cross_mask, cyto, save_folder):
    # Check if length of mask is equal to length of data frame
    checklen_mask = len(df_st['gene_symbol'][sig_13_mask]) + len(df_st['gene_symbol'][sig_24_mask]) +\
                    len(df_st['gene_symbol'][not_sig_mask]) + len(df_st['gene_symbol'][cross_mask])
    if (len(df_st['gene_symbol']) != checklen_mask) & (len(df_sc['gene_symbol']) != checklen_mask):
        raise ValueError

    # Check if mask dont include all points -> plot should be empty
    whole_mask = sig_13_mask | sig_24_mask | not_sig_mask | cross_mask
    fig, ax = plt.subplots(figsize=fig_size)
    ax.scatter(df_st[~whole_mask], df_sc[~whole_mask], c='k', s=dotsize, zorder=2, alpha=1, edgecolors='k',
               linewidth=0.2, label=r"Forgotten points")
    ax.legend()
    fig.savefig(os.path.join(save_folder, "_".join([cyto, "Unittest_Forgotten_points.pdf"])))
    plt.close()

    # Check if masks dont overlay
    whole_mask = sig_13_mask & sig_24_mask & not_sig_mask & cross_mask
    fig, ax = plt.subplots(figsize=fig_size)
    ax.scatter(df_st[whole_mask], df_sc[whole_mask], c='k', s=dotsize, zorder=2, alpha=1, edgecolors='k',
               linewidth=0.2, label=r"Intersection points")
    ax.legend()
    fig.savefig(os.path.join(save_folder, "_".join([cyto, "Unittest_Intersection_points.pdf"])))
    plt.close()

    # Check combinations of masks -> they should not overlapp --> if they do not overlap print 0
    m_allquadrants = sig_13_mask & sig_24_mask
    print(np.count_nonzero(m_allquadrants))
    m_13_notsig = sig_13_mask & not_sig_mask
    print(np.count_nonzero(m_13_notsig))
    m_13_cross = sig_13_mask & cross_mask
    print(np.count_nonzero(m_13_cross))
    m_24_notsig = sig_24_mask & not_sig_mask
    print(np.count_nonzero(m_24_notsig))
    m_24_cross = sig_24_mask & cross_mask
    print(np.count_nonzero(m_24_cross))
    m_notsig_cross = not_sig_mask & cross_mask
    print(np.count_nonzero(m_notsig_cross))

    # Plot all points -> all black points should be over-plotted by coloured ones
    fig, ax = plt.subplots(figsize=fig_size)
    # not significant: pval > 0.05
    ax.scatter(df_st[value], df_sc[value], c='k', s=dotsize, zorder=2, alpha=1, edgecolors='k', linewidth=0.2,
               label=r"All points")
    ax.scatter(df_st[value][cross_mask], df_sc[value][cross_mask], c='lightgrey', s=dotsize, zorder=2, alpha=1,
               edgecolors='lightgrey', linewidth=0.2, label=r"p-value > 0.05")
    # not significant and log2fc < threshold
    ax.scatter(df_st[value][not_sig_mask], df_sc[value][not_sig_mask], c='darkorange', s=dotsize, alpha=1,
               edgecolors='darkorange', linewidth=0.2, zorder=2, label=r"p-value < 0.05 and |log$_2$FC| < 1")

    # Add alpha to those which are complementary regulated
    ax.scatter(df_st[value][sig_24_mask], df_sc[value][sig_24_mask], c='darkred', s=dotsize,
               zorder=2, edgecolors='k', alpha=0.8, linewidth=0.2, label="p-value < 0.05 and |log$_2$FC| > 1")
    # Plot commonly significant genes: p-value < 0.05 and |log$_2$FC| > 1
    ax.scatter(df_st[value][sig_13_mask], df_sc[value][sig_13_mask], c='darkred', s=dotsize, zorder=2,
               edgecolors='k', alpha=1, linewidth=0.2, label=r"p-value < 0.05 and |log$_2$FC| > 1")
    ax.legend()
    fig.savefig(os.path.join(save_folder, "_".join([cyto, "Unittest_plot.pdf"])))
    plt.close()


def get_updowninbetween_masks(df_st, df_sc, significance_level, value, log2fc_cut=1., threshold=0.05):
    """Create masks for all scenarios (inside and outside of cross)

    Parameters
    ----------
    df_st : pandas.Dataframe
    df_sc : pandas.Dataframe
    significance_level: str
    pval or padj
    value : str
        name of p-value or FDR value to plot
    log2fc_cut : float
    threshold : float

    Returns
    -------

    """

    # 1. Significant p <= 0.05 (alpha darkorange)
    m_outsidecross_st = df_st[significance_level].values <= threshold
    m_outsidecross_sc = df_sc[significance_level].values <= threshold
    # -> all 4 quadrants: I, II, III, IV outside cross p < 0.05
    outsidecross_mask = np.logical_and(m_outsidecross_st, m_outsidecross_sc)

    # 2. Not Significant p > 0.05 (alpha grey)
    m_notpval_st = df_st[significance_level].values > threshold
    m_notpval_sc = df_sc[significance_level].values > threshold
    # -> all 4 quadrants in cross p > 0.05
    cross_mask = np.logical_or(m_notpval_st, m_notpval_sc)

    # 3. Significant p <= 0.05 log2fc >= 1 (darkred; balk circle)
    m_sig3_st = (df_st[value].values <= np.log10(threshold)) & (abs(df_st['log2fc'].values) >= log2fc_cut)
    m_sig3_sc = (df_sc[value].values <= np.log10(threshold)) & (abs(df_sc['log2fc'].values) >= log2fc_cut)
    m_sig1_st = (df_st[value].values >= -np.log10(threshold)) & (abs(df_st['log2fc'].values) >= log2fc_cut)
    m_sig1_sc = (df_sc[value].values >= -np.log10(threshold)) & (abs(df_sc['log2fc'].values) >= log2fc_cut)
    # -> quadrants I and III outside of cross
    sig_13_mask = (m_sig3_st & m_sig3_sc) | (m_sig1_st & m_sig1_sc)

    # 4. Complementary Significant (alpha darkred)
    m_sig2_st = (df_st[value].values >= -np.log10(threshold)) & (abs(df_st['log2fc'].values) >= log2fc_cut)
    m_sig2_sc = (df_sc[value].values <= np.log10(threshold)) & (abs(df_sc['log2fc'].values) >= log2fc_cut)
    m_sig4_st = (df_st[value].values <= np.log10(threshold)) & (abs(df_st['log2fc'].values) >= log2fc_cut)
    m_sig4_sc = (df_sc[value].values >= -np.log10(threshold)) & (abs(df_sc['log2fc'].values) >= log2fc_cut)
    # -> quadrants II and IV outside of cross
    sig_24_mask = (m_sig2_st & m_sig2_sc) | (m_sig4_st & m_sig4_sc)

    # 5. Outside cross without Significant: |log2FC| < 1 but p < 0.05
    not_sig_mask = outsidecross_mask & ~sig_13_mask & ~sig_24_mask

    index_hkg_st = np.where(df_st['gene_symbol'][:, np.newaxis] == np.array(_house_keeping_genes)[np.newaxis, :])[0]
    index_hkg_sc = np.where(df_sc['gene_symbol'][:, np.newaxis] == np.array(_house_keeping_genes)[np.newaxis, :])[0]

    return sig_13_mask, sig_24_mask, not_sig_mask, cross_mask, index_hkg_st, index_hkg_sc


def plotly_interactive_singedppvalues(df_st, df_sc, significance_level, value, save_folder, key, data_sets,
                                      log2fc_cut=1., threshold=0.05, zoom=True):
    """Plot interactive plot using plotly

    Parameters
    ----------
    df_st : pandas.Dataframe
        Spatial DGE Analysis results
    df_sc : pandas.Dataframe
        Single cell DGE Analysis results
    significance_level : str
        pval or padj
    value : str
    save_folder : str
        path to save folder
    key : str
    data_sets : 'list' ['str']
        contains names of data sets
    log2fc_cut : float
        threshold for effect size / log2FC
    threshold : float
        threshold for p-value and p-adjusted value
    zoom : bool
        if signature cytokine is removed or not

    Returns
    -------

    """
    # create masks for all scenarios
    sig_13_mask, sig_24_mask, not_sig_mask, cross_mask, index_hkg_st, index_hkg_sc = get_updowninbetween_masks(
        df_st=df_st, df_sc=df_sc, threshold=threshold, log2fc_cut=log2fc_cut,
        value=value, significance_level=significance_level)
    # Check masks
    check_masks(df_st, df_sc, value, sig_13_mask, sig_24_mask, not_sig_mask, cross_mask, cytokine, save_folder)

    if 'p' in value:
        xy_labels = r'signed log$_{10}$p-values'
        legend_label = 'p-value'
    else:
        xy_labels = r'signed log$_{10}$FDR-corrected p-value'
        legend_label = 'FDR'

    # plot
    fig = go.Figure()
    # not significant: pval > 0.05
    fig.add_trace(go.Scatter(x=df_st[value][cross_mask], y=df_sc[value][cross_mask],
                             mode='markers', marker=dict(size=dotsize, color='lightgrey', opacity=1),
                             text=df_st['gene_symbol'][cross_mask],  # hover text goes here
                             name=r"{} > 0.05".format(legend_label)))

    # not significant and log2fc < threshold
    fig.add_trace(go.Scatter(x=df_st[value][not_sig_mask], y=df_sc[value][not_sig_mask],
                             mode='markers', marker=dict(size=dotsize, color='darkorange', opacity=0.8),
                             text=df_st['gene_symbol'][not_sig_mask],  # hover text goes here
                             name=r"{} < 0.05 and |log$_2$FC| < 1".format(legend_label)))

    # plot significant genes
    fig.add_trace(go.Scatter(x=df_st[value][sig_13_mask], y=df_sc[value][sig_13_mask],
                             mode='markers', marker=dict(size=dotsize, color='darkred', opacity=1),
                             text=df_st['gene_symbol'][sig_13_mask],
                             name=r'{} <= 0.05 and |log$_2$FC| >= 1'.format(legend_label)))

    fig.add_trace(go.Scatter(x=df_st[value][sig_24_mask], y=df_sc[value][sig_24_mask],
                             mode='markers', marker=dict(size=dotsize, color='darkred', opacity=0.5),
                             text=df_st['gene_symbol'][sig_24_mask], showlegend=False))

    # Mark house keeping genes
    fig.add_trace(go.Scatter(x=df_st[value].iloc[index_hkg_st], y=df_sc[value].iloc[index_hkg_sc],
                             marker=dict(size=dotsize, color='darkgreen'), mode='markers',
                             text=df_st['gene_symbol'].iloc[index_hkg_st], name='Housekeeping genes'))

    fig.update_layout(title=xy_labels.upper() + " " + key, xaxis_title="{} \n {}".format(xy_labels, data_sets[0]),
                      yaxis_title="{} \n {}".format(xy_labels, data_sets[1]),
                      font=dict(family="Courier New, monospace", size=18, color="#7f7f7f"))

    # fig.write_image(os.path.join(save_folder, key + "_non_interactive.png"))
    if zoom:
        fig.write_html(os.path.join(save_folder, "_".join(["Zoom", key, "interactive.html"])))
    else:
        fig.write_html(os.path.join(save_folder, "_".join(["Full", key, "interactive.html"])))


def plot_signed_ppvalues(df_st, df_sc, significance_level, label_genes, save_folder, value, data_sets, sig_r,
                         threshold=0.05, adjust=False, zoom=True):
    """Plot signed and log10 transformed p-values

    Parameters
    ----------
    df_st : pandas.Dataframe
    df_sc : pandas.Dataframe
    significance_level : str
    label_genes : dict
    save_folder : str
    value : str
    data_sets : 'list' ['str']
    sig_r : pandas.Dataframe
    threshold : float
    adjust : bool
    zoom : bool

    Returns
    -------

    """

    # Get masks
    sig_13_mask, sig_24_mask, not_sig_mask, cross_mask, _, _ = get_updowninbetween_masks(
        df_st=df_st, df_sc=df_sc, threshold=threshold, log2fc_cut=1, value=value, significance_level=significance_level)
    # Check masks
    check_masks(df_st, df_sc, value, sig_13_mask, sig_24_mask, not_sig_mask, cross_mask, cytokine, save_folder)

    if 'p' in value:
        xy_labels = r'signed log$_{10}$p-values'
        legend_label = 'p-value'
    else:
        xy_labels = r'signed log$_{10}$FDR-corrected p-value'
        legend_label = 'FDR'

    fig, ax = plt.subplots(figsize=fig_size)
    # mark cross in signed p-value plot
    ax.axhline(y=0, c='grey', linestyle='dashed', lw=1, zorder=1)
    ax.axvline(x=0, c='grey', linestyle='dashed', lw=1, zorder=1)
    ax.axhline(y=np.log10(threshold), c='black', linestyle='dashed', lw=0.5, zorder=1)
    ax.axhline(y=-np.log10(threshold), c='black', linestyle='dashed', lw=0.5, zorder=1)
    ax.axvline(x=-np.log10(threshold), c='black', linestyle='dashed', lw=0.5, zorder=1)
    ax.axvline(x=np.log10(threshold), c='black', linestyle='dashed', lw=0.5, zorder=1)
    ax.set_xlabel("{} \n {}".format(xy_labels, data_sets[0]), fontsize=xy_fontsize)
    ax.set_ylabel("{} \n {}".format(xy_labels, data_sets[1]), fontsize=xy_fontsize)

    # plot data points
    # not significant: pval > 0.05
    ax_notsig = ax.scatter(df_st[value][cross_mask], df_sc[value][cross_mask], c='lightgrey', s=dotsize, zorder=2,
                           alpha=0.5, edgecolors='lightgrey', linewidth=0.2, label=r"{} > 0.05".format(legend_label))
    # not significant and log2fc < threshold
    ax_notsig_noteffectsize = ax.scatter(df_st[value][not_sig_mask], df_sc[value][not_sig_mask], c='darkorange',
                                         s=dotsize, alpha=0.5, edgecolors='darkorange', linewidth=0.2, zorder=2,
                                         label=r"{} < 0.05 and |log$_2$FC| < 1".format(legend_label))

    # Add alpha to those which are complementary regulated
    ax.scatter(df_st[value][sig_24_mask], df_sc[value][sig_24_mask], c='darkred', s=dotsize, zorder=2,
               edgecolors='k', alpha=0.8, linewidth=0.2,
               label=r"{} <= 0.05 and |log$_2$FC| >= 1".format(legend_label))
    # Plot commonly significant genes: p-value < 0.05 and |log$_2$FC| > 1
    ax_sig = ax.scatter(df_st[value][sig_13_mask], df_sc[value][sig_13_mask], c='darkred', s=dotsize, zorder=2,
                        edgecolors='k', alpha=1, linewidth=0.2,
                        label=r"{} <= 0.05 and |log$_2$FC| >= 1".format(legend_label))

    # Add driver and responder gene annotations
    color_genes = get_colors_genes(df_st=df_st, df_sc=df_sc, label_genes=label_genes, value=value)
    texts = []
    for ind_label, (x, y) in enumerate(zip(color_genes['x'], color_genes['y'])):
        texts.append(ax.text(x, y, color_genes['gene_symbol'].values[ind_label], size=text_fontsize,
                             color=color_genes['color'].values[ind_label], zorder=3))

    if adjust:
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=1, zorder=3), ax=ax, lim=1000,
                    expand_text=(2, 2))

    # Add Threshold
    log10_cut = -np.log10(threshold)
    if cytokine == 'IL17A':
        # Full: 83; Zoomed: 39
        if zoom:
            ax.text(33, log10_cut + 0.3, "5% FDR", size=10, color='k', zorder=3)
            ax.text(log10_cut, np.amin(df_sc[value].values), "5% FDR", size=10, color='k', zorder=3)
            ax.set_ylim([-15, 50])
            ax.set_xlim([-28, 38])
            # Add correlation and correlation p-value
            ax.text(np.amax(df_st[value].values) / 2 - np.amax(df_st[value].values) / 10, 51,
                    'r = {:.2f}; p = {:.2e}'.format(sig_r['r'].values[0], sig_r['p-val'].values[0]),
                    fontstyle='italic', fontsize=text_fontsize)
        else:
            ax.text(85, log10_cut + 0.2, "5% FDR", size=8, color='k', zorder=3)
            ax.text(log10_cut + 0.2, np.amin(df_sc[value].values) - 0.5, "5% FDR", size=8, color='k', zorder=3)
            # Add correlation and correlation p-value
            ax.text(np.amax(df_st[value].values) / 2 - np.amax(df_st[value].values) / 10,
                    np.amax(df_sc[value].values) + np.amax(df_sc[value].values) / 8,
                    'r = {:.2f}; p = {:.2e}'.format(sig_r['r'].values[0], sig_r['p-val'].values[0]),
                    fontstyle='italic', fontsize=text_fontsize)
    else:
        if zoom:
            ax.set_ylim([-18, 25])
            ax.set_xlim([-18, 48])
            ax.text(43, log10_cut + 0.1, "5% FDR", size=8, color='k', zorder=3)
            ax.text(log10_cut + 0.2, np.amin(df_sc[value].values), "5% FDR", size=8, color='k', zorder=3)
        else:
            ax.text(np.amax(df_st[value].values) - 4, log10_cut + 0.5, "5% FDR", size=8, color='k', zorder=3)
            ax.text(log10_cut + 0.2, np.amin(df_sc[value].values) - 0.5, "5% FDR", size=8, color='k', zorder=3)
        # Add correlation and correlation p-value
        ax.text(np.amax(df_st[value].values) / 2 - np.amax(df_st[value].values) / 10,
                np.amax(df_sc[value].values) + np.amax(df_sc[value].values) / 8,
                'r = {:.2f}; p = {:.2e}'.format(sig_r['r'].values[0], sig_r['p-val'].values[0]),
                fontstyle='italic', fontsize=text_fontsize)

    # Add legend, Labels Driver and Responder genes
    # create blank rectangle
    driver_rect = Rectangle((0, 0), 1, 1, fc="w", fill=True, edgecolor='w', linewidth=0)
    resp_rect = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='w', linewidth=0)
    leg = ax.legend([ax_sig, ax_notsig_noteffectsize, ax_notsig, driver_rect, resp_rect],
                    (r"{} $\leq$ 0.05 and |log$_2$FC| $\geq$ 1".format(legend_label),
                     r"{} $<$ 0.05 and |log$_2$FC| $<$ 1".format(legend_label),
                     "{} $>$ 0.05".format(legend_label), "Leukocyte genes", "{} responder genes".format(cytokine)))
    leg_color = ['k', 'k', 'k', 'purple', 'mediumblue']
    for ind, text in enumerate(leg.get_texts()):
        plt.setp(text, color=leg_color[ind])

    # remove borders from plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)

    plt.tight_layout()

    if zoom:
        fig.savefig(os.path.join(save_folder, "_".join(["_".join(['Zoomed', cytokine, value, "plot"]),
                                                        dge_method, design_function, '.pdf'])))
    else:
        fig.savefig(os.path.join(save_folder, "_".join(["_".join(['Full', cytokine, value, "plot"]),
                                                        dge_method, design_function, '.pdf'])))

    plt.close()


def get_colors_genes(df_st, df_sc, label_genes, value):
    """Get colors for driver and responder genes

    This function creates a data frame for driver and responder genes containing the signed p-values as x and y
    and add the colors purple and blue for driver and responder genes, respectively.

    Parameters
    ----------
    df_st : pandas.Dataframe
        Data frame containing the results of the dge analysis
    df_sc : pandas.Dataframe
        Data frame containing the results of the dge analysis
    label_genes : dict
        contains the gene symbols for driver and responder genes
    value : str
        Column name of value in df

    Returns
    -------
    pandas.Dataframe
        Contains a genes' position (x=ST p-value, y=SC p-value)

    Examples
    --------
    """
    label_info = pd.DataFrame()

    cytokines = ['IL17A', 'IFNG', 'IL13']

    driver_group = label_genes['Driver_genes']
    responder_group = label_genes['Responder_genes']
    merged_genes = driver_group.copy()
    merged_genes.extend(responder_group)
    merged_genes.extend(cytokines)

    # Get index of driver and responder genes from ST and SC dge analysis results
    ind_st_genes = np.where(df_st['gene_symbol'][np.newaxis, :] == np.array(merged_genes)[:, np.newaxis])[1]
    ind_sc_genes = np.where(df_sc['gene_symbol'][np.newaxis, :] == np.array(merged_genes)[:, np.newaxis])[1]
    label_st_df = df_st.iloc[ind_st_genes]
    label_sc_df = df_sc.iloc[ind_sc_genes]

    # Get singed p-values from ST and SC data set as x and y position
    label_info['x'], label_info['y'] = label_st_df[value].values, label_sc_df[value].values
    label_info['gene_symbol'] = label_st_df['gene_symbol'].values
    # Assign color for driver and responder genes
    label_info['color'] = ['mediumblue'] * len(label_info.index)
    ind_drivers = np.where(label_info['gene_symbol'][np.newaxis, :] == np.array(driver_group)[:, np.newaxis])[1]
    label_info['color'].values[ind_drivers] = ['purple'] * len(ind_drivers)

    # Add color for cytokines
    ind_cytos = np.where(label_info['gene_symbol'][np.newaxis, :] == np.array(cytokines)[:, np.newaxis])[1]
    if cytokine == 'IL17A':
        label_info['color'].values[ind_cytos] = ['#ff7f00'] * len(ind_cytos)
    else:
        label_info['color'].values[ind_cytos] = ['#377eb8'] * len(ind_cytos)

    return label_info


def get_correlation(df_st, df_sc, significance_level, value, pval_cut=0.05, log2fc_cut=1., method='pearson'):
    """Calculate Pearson Correlation

    Parameters
    ----------
    df_st : pandas.Dataframe
        Data frame containing the results of the dge analysis
    df_sc : pandas.Dataframe
        Data frame containing the results of the dge analysis
    significance_level : str
        pval or padj
    value : str
        name of signed value
    pval_cut : float
        p-value cut parameter
    log2fc_cut : float
        effect size cut parameter
    method : str
        Correlation method to use

    Returns
    -------

    """
    # 1. get up-regulated genes in cyto+ group in both data sets
    m_sig_st = (df_st[significance_level].values <= pval_cut) & (abs(df_st['log2fc'].values) >= log2fc_cut)
    m_sig_sc = (df_sc[significance_level].values <= pval_cut) & (abs(df_sc['log2fc'].values) >= log2fc_cut)
    mask_siggenes = np.logical_and(m_sig_st, m_sig_sc)

    # 2. Calculate Pearson Correlation and p-value
    sig_r = pingouin.corr(x=df_st[value][mask_siggenes], y=df_sc[value][mask_siggenes],
                          method=method)
    print(sig_r['p-val'].values[0])
    print(sig_r['r'])

    return sig_r


def get_index(df_st, df_sc, gene_set):
    """Get index of genes in gene_set

    Parameters
    ----------
    df_st : pandas.Dataframe
        Data frame containing the results of the dge analysis
    df_sc : pandas.Dataframe
        Data frame containing the results of the dge analysis
    gene_set : set
        genes of interest

    Returns
    -------

    """
    st_index_genes = np.where(np.array(df_st['gene_symbol'])[:, np.newaxis] ==
                              np.array(list(gene_set))[np.newaxis, :])[0]
    sc_index_genes = np.where(np.array(df_sc['gene_symbol'])[:, np.newaxis] ==
                              np.array(list(gene_set))[np.newaxis, :])[0]

    return st_index_genes, sc_index_genes


def load_data(path):
    """Load data by given input path

    Parameters
    ----------
    path : str
        Path to input data

    Returns
    -------
    pandas.Dataframe
    """
    df_data = pd.read_csv(path, error_bad_lines=False)

    # Remove column Unnamed: 0
    df_data = df_data.drop(['Unnamed: 0'], axis=1)

    # Check if column names and row names are unique
    print("Unique Genes data:", df_data['gene_symbol'].is_unique)
    print("Unique Columns data:", df_data.columns.is_unique)

    # remove duplicated rows
    df_data = df_data.loc[~df_data['gene_symbol'].duplicated(), :]

    return df_data


def data_preparation(path_st_data, path_sc_data):
    """

    Parameters
    ----------
    path_st_data : str
        path to DGE Analysis result of the Spatial transcriptomic data set
    path_sc_data : str
        path to DGE Analysis result of the Single-cell transcriptomic data set

    Returns
    -------

    """
    # load files
    df_st_data = load_data(path=path_st_data)
    df_sc_data = load_data(path=path_sc_data)

    # get common genes
    genes = set(df_st_data['gene_symbol'].values) & set(df_sc_data['gene_symbol'].values)
    # get index of common genes
    index_st, index_sc = get_index(df_st=df_st_data, df_sc=df_sc_data, gene_set=genes)
    # subset df's to common genes
    df_st_data = df_st_data.iloc[index_st]
    df_sc_data = df_sc_data.iloc[index_sc]

    # sort rows of ST data by SC data frame
    df_st_data = df_st_data.set_index('gene_symbol')
    df_st_data = df_st_data.reindex(index=df_sc_data['gene_symbol'])
    df_st_data = df_st_data.reset_index()

    return df_st_data, df_sc_data


def main(path_st_data, path_sc_data, save_folder, log2fc_cut=1.0, pval_cut=0.05, zoom=True):
    """Plot signed log10-transformed p-values of ST vs. scRNA-seq DGE Analysis results

    Parameters
    ----------
    path_st_data : str
        path to DGE Analysis result of the Spatial transcriptomics data set
    path_sc_data : str
        path to DGE Analysis result of the Single-cell transcriptomics data set
    save_folder : str
        path to result storing folder
    log2fc_cut : float
    pval_cut : float
    zoom : bool

    Returns
    -------

    """
    print("1. Get common genes between ST and SC DGE Analysis results")
    df_st, df_sc = data_preparation(path_st_data=path_st_data, path_sc_data=path_sc_data)

    # Drop cytokine since it has the highest p-value and was the selection criteria
    if zoom:
        df_st = df_st.drop([0])
        df_sc = df_sc.drop([0])

    # Switch signs - new since 26.09.2021
    # df_st['log2fc'] = -df_st['log2fc']
    # df_sc['log2fc'] = -df_sc['log2fc']

    # Get signed log10(p-values)
    df_st['signed_FDR'] = np.sign(df_st['log2fc']) * np.log10(df_st['padj'])
    df_sc['signed_FDR'] = np.sign(df_sc['log2fc']) * np.log10(df_sc['padj'])
    print("2. Get driver and responder genes")
    genes = gene_lists.highlight_genes()

    # Enrichment test of significant value + their Correlation
    # -> use Spearman Correlation because its not so sensitive towards outliers
    sig_corr = get_correlation(df_st=df_st, df_sc=df_sc, significance_level='padj', value='signed_FDR',
                               pval_cut=pval_cut, log2fc_cut=log2fc_cut, method='spearman')

    print("3. Plot signed pp-value plot")
    plot_signed_ppvalues(df_st=df_st, df_sc=df_sc, label_genes=genes[cytokine], save_folder=save_folder,
                         significance_level='padj', value='signed_FDR', data_sets=['ST dataset', 'scRNA-seq dataset'],
                         threshold=pval_cut, adjust=True, sig_r=sig_corr, zoom=zoom)

    # Interactive signed p-value plot
    plotly_interactive_singedppvalues(
        df_st=df_st, df_sc=df_sc, value='signed_FDR', save_folder=save_folder, key=cytokine, significance_level='padj',
        data_sets=['ST dataset', 'scRNA-seq dataset'], log2fc_cut=log2fc_cut, threshold=pval_cut, zoom=zoom)


if __name__ == '__main__':
    # 0. Initialisation
    cytokine = 'IL17A'
    design_function = 'cdr_patient_annotation_cyto'
    dataset = "Whole_T_cell_matrix"
    dge_method = 'glmGamPoi'
    log2fc_threshold = 1
    pval_threshold = 0.05
    zoom_in = True

    # 1. get date
    today = date.today()

    # 2. create saving folder in current project path
    savepath = os.path.join("..", "..", "..", "output", "Figure_3E", str(today))
    os.makedirs(savepath, exist_ok=True)

    # 3. Load DGE Analysis results
    input_path = os.path.join("..", "..", "..", "input", "dge_analysis")

    # load df's using pandas
    comparisons = "_".join([cytokine, 'vs_Others'])
    st_path = os.path.join(
        input_path, "".join(['2021-02-01_spatial__cdr_patient_annotation_cyto', os.path.sep, 'spatial',
                             os.path.sep, dataset, "__", design_function,
                             os.path.sep, cytokine, os.path.sep, comparisons, os.path.sep,
                             "spatial_", comparisons, "_glmGamPoi_DGE_all_genes.csv"]))

    sc_path = os.path.join(
        input_path, "".join(['2021-02-01_single_cell__cdr_annotation_cyto', os.path.sep, 'single_cell', os.path.sep,
                             dataset, "__", "cdr_annotation_cyto", os.path.sep, cytokine,
                             os.path.sep, comparisons, os.path.sep, "single_cell_", comparisons,
                             "_glmGamPoi_DGE_all_genes.csv"]))

    main(path_st_data=st_path, path_sc_data=sc_path, save_folder=savepath,
         log2fc_cut=log2fc_threshold, pval_cut=pval_threshold, zoom=zoom_in)
