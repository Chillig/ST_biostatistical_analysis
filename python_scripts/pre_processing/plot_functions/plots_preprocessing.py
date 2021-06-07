import numpy as np
import scanpy as sc
import os
import pandas as pd

import seaborn as sb
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

import python_scripts.pre_processing.plot_functions.helper_functions as hf

# Initializing scanpy
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
try:
    sc.logging.print_versions()
except AssertionError:
    sc.logging.print_versions()
# sc.settings.set_figure_params(dpi=80)
matplotlib.get_backend()

# containment radius
_1sigma = 0.682 * 100
_1_5sigma = 0.866 * 100
_2sigma = 0.954 * 100
_3sigma = 0.997 * 100
_simga_string = "\u03C3"

fig_size, title_fontsize, subtitle_fontsize, xy_fontsize, xy_ticks, legend_fontsize, text_fontsize, fileformat = \
    hf.plot_properties()


def _plot_umap(amatrix, ax=None, color='sample', use_raw=False):
    """Basic UMAP plot

    Parameters
    ----------
    amatrix : annData
    ax : matplotlib.pyplot.axis
    color : str
    use_raw : bool
        if raw annData is provided

    Returns
    -------

    """
    img_ax = sc.pl.umap(amatrix, color=color, use_raw=use_raw, ax=ax, wspace=0.4, show=False)
    return img_ax


def plot_qc_metrics(amatrix, save_folder, sample_name, counts_threshold=10000, raw=False):
    """Sample quality plots

    Parameters
    ----------
    amatrix : annData
    save_folder : str
    sample_name : str
    counts_threshold : int
    raw : bool

    Returns
    -------

    """

    # Sample quality plots
    fig, axs = plt.subplots(nrows=1, ncols=1, facecolor='w', edgecolor='k', figsize=(18, 12))
    fig.subplots_adjust(hspace=1, wspace=0.3)
    # axs = axs.ravel()
    # 1. UMI counts
    sc.pl.violin(amatrix, 'n_counts', groupby='sample', size=2, log=True, cut=0, rotation=90, ax=axs, show=False)
    hf.set_axes_label(x_label='Sample', y_label='UMI counts', axes=axs)
    fig.suptitle("Sample quality plots", fontsize=16)
    fig.tight_layout(pad=1, w_pad=1, h_pad=1.0)
    fig.subplots_adjust(top=0.88)

    if raw:
        plt.savefig(os.path.join(save_folder, "_".join([sample_name, "raw_UMI_counts_Sample_QC.png"])))
    else:
        plt.savefig(os.path.join(save_folder, "_".join([sample_name, "UMI_counts_Sample_QC.png"])))

    fig, axs = plt.subplots(nrows=1, ncols=1, facecolor='w', edgecolor='k', figsize=(18, 12))
    fig.subplots_adjust(hspace=1, wspace=0.3)
    # 2. MT-fraction
    sc.pl.violin(amatrix, 'mt_frac', groupby='sample', ax=axs, rotation=90, show=False)
    hf.set_axes_label(x_label='Sample', y_label='Mitochondrial fraction [%]', axes=axs)
    fig.suptitle("Sample quality plots", fontsize=16)
    fig.tight_layout(pad=1, w_pad=1, h_pad=1.0)
    fig.subplots_adjust(top=0.88)

    if raw:
        plt.savefig(os.path.join(save_folder, "_".join([sample_name, "raw_MT-fraction_Sample_QC.png"])))
    else:
        plt.savefig(os.path.join(save_folder, "_".join([sample_name, "MT-fraction_Sample_QC.png"])))

    # Data quality summary plots
    fig1, axs1 = plt.subplots(nrows=1, ncols=2, facecolor='w', edgecolor='k', figsize=(12, 6))
    fig1.subplots_adjust(hspace=1, wspace=0.5)
    axs1 = axs1.ravel()
    sc.pl.scatter(amatrix, 'n_counts', 'n_genes', color='mt_frac', ax=axs1[0], right_margin=1.5, show=False)
    hf.set_axes_label(x_label='No. UMI counts', y_label='No. genes', axes=axs1[0])
    hf.set_title(title='Mitochondrial fraction', axes=axs1[0])
    sc.pl.scatter(amatrix[amatrix.obs['n_counts'] < counts_threshold], 'n_counts', 'n_genes', color='mt_frac',
                  ax=axs1[1], right_margin=0.12, show=False)
    hf.set_axes_label(x_label='No. UMI counts', y_label='No. genes', axes=axs1[1])
    hf.set_title(sup_title="Data quality summary plots", figure=fig1)
    hf.set_title(title='Mitochondrial fraction without outliers', axes=axs1[1])
    fig1.subplots_adjust(top=0.88)

    if raw:
        plt.savefig(os.path.join(save_folder, "_".join([sample_name, "raw_Data_QC.png"])))
    else:
        plt.savefig(os.path.join(save_folder, "_".join([sample_name, "Data_QC.png"])))
    plt.close()


def plot_distribution(amatrix, save_folder, sample_name, lower_filter_counts=400, upper_filter_counts=4000,
                      upper_filter_genes=500, log_scale=None, bins=np.ones(6) * 60, raw=False):
    """Plot distribution of counts
        Threshold decision based on the distribution of counts
        Distribution should be shaped like a t-distribution or gaussian distribution to apply later statistical tests

    Parameters
    ----------
    amatrix : annData
    save_folder : str
    sample_name : str
    lower_filter_counts : int
    upper_filter_counts : int
    upper_filter_genes : int
    log_scale : None, bool
    bins : numpy.array
    raw : bool

    Returns
    -------

    """
    # Plots for threshold determination of n_counts
    fig = plt.figure(facecolor='w', edgecolor='k', figsize=(10, 6))
    fig.subplots_adjust(bottom=0.025, left=0.025, top=0.975, right=0.975)
    axs = [(2, 2, (1, 2)), (2, 2, 3), (2, 2, 4)]
    if log_scale:
        for nrows, ncols, plot_number in axs:
            sub = fig.add_subplot(nrows, ncols, plot_number)
            sub.grid(which='both', linestyle='--', zorder=0)
            if plot_number == (1, 2):
                p3, bins0, _ = plt.hist(amatrix.obs['n_counts'], bins=bins[0], color='blue', edgecolor='k',
                                        linewidth=0.5, zorder=3)
                # calculate sigma containment radius
                hf.place_cut_txtbox(lower_filter_counts, ax=sub, y_pos_txtbox=80)
                hf.place_cut_txtbox(upper_filter_counts, ax=sub, y_pos_txtbox=80)
                binwidth = abs(np.amax(bins0) - np.amin(bins0)) / bins[0]
                hf.set_title(title='UMI count distribution with binwidth = {}'.format(int(binwidth)), axes=sub)
            elif plot_number == 3:
                p4, bins1, _ = plt.hist(amatrix.obs['n_counts'][amatrix.obs['n_counts'] < lower_filter_counts],
                                        bins=bins[1], color='blue', edgecolor='k', linewidth=0.5, zorder=3)
                binwidth = abs(np.amax(bins1) - np.amin(bins1)) / bins[1]
                hf.set_title(title='No. counts {} with binwidth = {}'.format(lower_filter_counts, int(binwidth)),
                             axes=sub)
            elif plot_number == 4:
                p5, bins2, _ = plt.hist(amatrix.obs['n_counts'][amatrix.obs['n_counts'] > upper_filter_counts],
                                        bins=bins[2], color='blue', edgecolor='k', linewidth=0.5, zorder=3)
                binwidth = abs(np.amax(bins2) - np.amin(bins2)) / bins[2]
                hf.set_title(title='Upper {} No. counts with binwidth = {}'.format(upper_filter_counts, int(binwidth)),
                             axes=sub)

            hf.set_axes_label(x_label='UMI counts', y_label='counts', axes=sub)
            sub.set_ylim([0.1, np.amax(p3) + np.amax(p3) * 0.05])
            sub.set_yscale('log', nonposy='clip')
    else:
        for nrows, ncols, plot_number in axs:
            sub = fig.add_subplot(nrows, ncols, plot_number)
            if plot_number == (1, 2):
                sb.distplot(amatrix.obs['n_counts'], kde=False, bins=bins[0], ax=sub)
            elif plot_number == 3:
                sb.distplot(amatrix.obs['n_counts'][amatrix.obs['n_counts'] < lower_filter_counts],
                            kde=False, bins=bins[1], ax=sub)
            elif plot_number == 4:
                sb.distplot(amatrix.obs['n_counts'][amatrix.obs['n_counts'] > upper_filter_counts],
                            kde=False, bins=bins[2], ax=sub)
            hf.set_axes_label(x_label='UMI counts', y_label='counts', axes=sub)
    hf.set_title(sup_title="Decision of UMI count threshold", figure=fig)
    fig.tight_layout(pad=1.5, w_pad=1.5, h_pad=1.5)
    fig.subplots_adjust(top=0.88)

    if raw:
        plt.savefig(os.path.join(save_folder, "_".join([sample_name, "raw_Distribution_counts.png"])))
    else:
        plt.savefig(os.path.join(save_folder, "_".join([sample_name, "Distribution_counts.png"])))

    # Plots for threshold determination of n_genes
    fig_genes = plt.figure(facecolor='w', edgecolor='k', figsize=(10, 6))
    fig_genes.subplots_adjust(bottom=0.025, left=0.025, top=0.975, right=0.975)
    axs = [(1, 2, 1), (1, 2, 2)]
    if log_scale:
        for nrows, ncols, plot_number in axs:
            sub_genes = fig_genes.add_subplot(nrows, ncols, plot_number)
            sub_genes.grid(which='both', linestyle='--', zorder=0)
            if plot_number == 1:
                p6, bins3, _ = plt.hist(amatrix.obs['n_genes'].values, bins=bins[3], color='blue', edgecolor='k',
                                        linewidth=0.5, zorder=3)
                # place threshold box
                hf.place_cut_txtbox(upper_filter_genes, ax=sub_genes)
                binwidth = abs(np.amax(bins3) - np.amin(bins3)) / bins[3]
                hf.set_title(title='Gene histogram with binwidth = {}'.format(int(binwidth)), axes=sub_genes)
            elif plot_number == 2:
                p7, bins4, _ = plt.hist(amatrix.obs['n_genes'][amatrix.obs['n_genes'].values < upper_filter_genes],
                                        bins=bins[4], color='blue', edgecolor='k', linewidth=0.5, zorder=3)
                binwidth = abs(np.amax(bins4) - np.amin(bins4)) / bins[4]
                hf.set_title(title='{}% of all genes with binwidth = {}' .format(upper_filter_genes, int(binwidth)),
                             axes=sub_genes)
            hf.set_axes_label(x_label='No. genes', y_label='counts', axes=sub_genes)
            sub_genes.set_ylim([0.1, np.amax(amatrix.obs['n_genes'].values)])
            sub_genes.set_yscale('log')
    else:
        for nrows, ncols, plot_number in axs:
            sub_genes = fig_genes.add_subplot(nrows, ncols, plot_number)
            if plot_number == 1:
                sb.distplot(amatrix.obs['n_genes'], kde=False, bins=bins[3], ax=sub_genes)
            elif plot_number == 2:
                sb.distplot(amatrix.obs['n_genes'][amatrix.obs['n_genes'] < upper_filter_genes], kde=False,
                            bins=bins[4], ax=sub_genes)
            hf.set_axes_label(x_label='No. genes', y_label='counts', axes=sub_genes)

    hf.set_title(sup_title="Decision of gene threshold", figure=fig_genes)
    fig_genes.tight_layout(pad=1.5, w_pad=1.5, h_pad=1.5)
    fig_genes.subplots_adjust(top=0.88)

    if raw:
        plt.savefig(os.path.join(save_folder, "_".join([sample_name, "raw_Distribution_genes.png"])))
    else:
        plt.savefig(os.path.join(save_folder, "_".join([sample_name, "Distribution_genes.png"])))
    plt.close()


def plot_sc_dist(adata, save_folder, bins=None, title=None, ax_xlabel=None, ax_ylabel=None, observable=None,
                 border_value=None, kde=False, n_rows=1, n_cols=1, color="blue", raw=False):
    """
    Plots the distribution of a anndata matrix
    Can be used for thresholding decision: counts and genes

    :param adata: [annData] or [Series, 1d-array, or list]
        Observed data. If this is a Series object with a ``name`` attribute,
        the name will be used to label the data axis.
    :param save_folder:
    :param observable: int, optional
        investigate observable of anndata matrix.
    :param border_value: array, optional
        numpy array: [0] = border value and [1] = sign to determine chosen start/end of interval.
    :param kde: bool, optional
        Whether to plot a gaussian kernel density estimate.
    :param bins: argument for matplotlib hist(), or None, optional
        Specification of hist bins. If unspecified, as reference rule is used
        that tries to find a useful default.
    :param n_rows: int, optional
        number of rows in figure (indirectly determines number of plots in figure).
    :param n_cols: int, optional
        determine the number of columns in figure (indirectly determines number of plots in figure).
    :param title: string, optional
        title of plot.
    :param ax_xlabel: string, optional
        lable of x-axis.
    :param ax_ylabel: string, optional
        lable of y-axis.
    :param color: by default blue
    :param raw:
    :return: ---
    """

    if not bool(ax_ylabel) and not bool(ax_xlabel):
        ax_label = True
    else:
        ax_label = False

    fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=fig_size)

    if bool(observable):
        if bool(border_value):
            if border_value[1] == "greater_than":
                sb.distplot(adata.obs[observable][adata.obs[observable] > border_value[0]], kde=kde, bins=bins,
                            axlabel=ax_label, ax=ax, color=color)
            elif border_value[1] == "greater_equal_than":
                sb.distplot(adata.obs[observable][adata.obs[observable] >= border_value[0]], kde=kde, bins=bins,
                            axlabel=ax_label, ax=ax)
            elif border_value[1] == "equal":
                sb.distplot(adata.obs[observable][adata.obs[observable] == border_value[0]], kde=kde, bins=bins,
                            axlabel=ax_label, ax=ax)
            elif border_value[1] == "smaller_equal_than":
                sb.distplot(adata.obs[observable][adata.obs[observable] <= border_value[0]], kde=kde, bins=bins,
                            axlabel=ax_label, ax=ax)
            elif border_value[1] == "smaller_than":
                sb.distplot(adata.obs[observable][adata.obs[observable] < border_value[0]], kde=kde, bins=bins,
                            axlabel=ax_label, ax=ax)
        else:
            sb.distplot(adata.obs[observable], kde=kde, bins=bins, axlabel=ax_label, ax=ax)
    else:
        sb.distplot(adata, bins=bins, kde=kde, axlabel=ax_label, ax=ax)

    if bool(title):
        ax.set_title(title)
    if bool(ax_xlabel):
        ax.set_xlabel(ax_xlabel)
    if bool(ax_ylabel):
        ax.set_ylabel(ax_ylabel)

    if raw:
        plt.savefig(os.path.join(save_folder, "raw_dist_threshold.png"))
    else:
        plt.savefig(os.path.join(save_folder, "dist_threshold.png"))
    plt.close()


def plot_sc_scatter(adata, save_folder, obs, title=None, ax_xlabel=None, ax_ylabel=None, n_rows=1, n_cols=1, raw=False):
    """Uses Scanpy's scatter plot and offers the possibility to create plots in one figure

    Parameters
    ----------
    adata : annData
        Observed data. If this is a Series object with a ``name`` attribute,
        the name will be used to label the data axis.
    save_folder : str
    obs : list, str
        investigate observable of matrix.
    title : None, str
        title of plot
    ax_xlabel : None, str, list
        x-axis label
    ax_ylabel : None, str, list
        y-axis label
    n_rows : int
        number of rows in figure
    n_cols : int
        number of columns in figure
    raw : bool
        if raw annData is provided

    Returns
    -------

    """

    fig, axs = plt.subplots(nrows=n_rows, ncols=n_cols, facecolor='w', edgecolor='k', figsize=fig_size)
    fig.subplots_adjust(hspace=1, wspace=0.5)
    plt.title(title)
    if n_cols > 1 or n_rows > 1:
        axs = axs.ravel()
        for obs_counter in range(len(obs)-1):
            sc.pl.scatter(adata, obs[0], obs[obs_counter + 1], ax=axs[obs_counter])
            hf.set_axes_label(x_label=ax_xlabel, y_label=ax_ylabel[obs_counter], axes=axs[obs_counter])
    else:
        sc.pl.scatter(adata, obs[0], obs[1], ax=axs)
        hf.set_axes_label(x_label=ax_xlabel, y_label=ax_ylabel, axes=axs)
    plt.show()

    if raw:
        plt.savefig(os.path.join(save_folder, "raw_sizefactor.png"))
    else:
        plt.savefig(os.path.join(save_folder, "sizefactor.png"))
    plt.close()


def plot_visualization_results(adata, save_folder, batch_key, raw=False):
    """Visualise data in PC, UMAP and tsne  space

    Parameters
    ----------
    adata : annData
    save_folder : str
        path saving directory
    batch_key : str
        if batch correction was applied or not
    raw : bool
        if raw adata or filtered adata object

    Returns
    -------

    """
    if batch_key == "batch_corrected":
        fig_2d_visu = plt.figure(facecolor='w', edgecolor='k', figsize=(10, 6))
        fig_2d_visu.subplots_adjust(bottom=0.025, left=0.025, top=0.975, right=0.975)
        # axs = [(3, 2, 1), (3, 2, 2), (3, 2, 3), (3, 2, 4), (3, 2, 5), (3, 2, 6)]
        sub = fig_2d_visu.add_subplot(1, 2, 1)
        sc.pl.pca_scatter(adata, color='n_counts', ax=sub)
        hf.set_axes_label(x_label='PC_1', y_label='PC_2', axes=sub)

        sub3 = fig_2d_visu.add_subplot(1, 2, 2)
        _plot_umap(adata, ax=sub3, color='n_counts')
        hf.set_axes_label(x_label='UMAP_1', y_label='UMAP_2', axes=sub3)

    else:
        fig_2d_visu = plt.figure(facecolor='w', edgecolor='k', figsize=(10, 6))
        fig_2d_visu.subplots_adjust(bottom=0.025, left=0.025, top=0.975, right=0.975)
        # axs = [(3, 2, 1), (3, 2, 2), (3, 2, 3), (3, 2, 4), (3, 2, 5), (3, 2, 6)]
        sub = fig_2d_visu.add_subplot(2, 3, 1)
        sc.pl.pca_scatter(adata, color='n_counts', ax=sub, show=False)
        hf.set_axes_label(x_label='PC_1', y_label='PC_2', axes=sub)

        sub2 = fig_2d_visu.add_subplot(2, 3, 2)
        sc.pl.tsne(adata, color='n_counts', ax=sub2, show=False)
        hf.set_axes_label(x_label='tSNE_1', y_label='t-SNE_2', axes=sub2)

        sub3 = fig_2d_visu.add_subplot(2, 3, 3)
        _plot_umap(adata, ax=sub3, color='n_counts')
        hf.set_axes_label(x_label='UMAP_1', y_label='UMAP_2', axes=sub3)

        sub4 = fig_2d_visu.add_subplot(2, 3, 4)
        sc.pl.diffmap(adata, color='n_counts', ax=sub4, components=[1, 2], show=False)
        hf.set_axes_label(x_label='DC_1', y_label='DC_2', axes=sub4)

        sub5 = fig_2d_visu.add_subplot(2, 3, 5)
        sc.pl.diffmap(adata, color='n_counts', ax=sub5, components=[1, 3], show=False)
        hf.set_axes_label(x_label='DC_1', y_label='DC_2', axes=sub5)

        sub6 = fig_2d_visu.add_subplot(2, 3, 6)
        sc.pl.draw_graph(adata, color='n_counts', ax=sub6, show=False)
        hf.set_axes_label(x_label='FA_1', y_label='FA_2', axes=sub6)

    hf.set_title(sup_title="2D Data Visualization", figure=fig_2d_visu)
    fig_2d_visu.tight_layout(pad=1.5, w_pad=1.5, h_pad=1.5)
    fig_2d_visu.subplots_adjust(top=0.88)

    plt.show()
    if raw:
        plt.savefig(os.path.join(save_folder, "raw_2D_Visualisation_{}.png".format(batch_key)))
    else:
        plt.savefig(os.path.join(save_folder, "2D_Visualisation_{}.png".format(batch_key)))
    plt.close()


def plot_cell_cycle_cluster(adata, save_folder, raw=False):
    """Plot cell cycle phases

    Parameters
    ----------
    adata : annData
    save_folder : str
    raw : bool

    Returns
    -------

    """

    fig_2d_visu = plt.figure(facecolor='w', edgecolor='k', figsize=(10, 6))
    fig_2d_visu.subplots_adjust(bottom=0.025, left=0.025, top=0.975, right=0.975)
    # axs = [(3, 2, 1), (3, 2, 2), (3, 2, 3), (3, 2, 4), (3, 2, 5), (3, 2, 6)]
    sub = fig_2d_visu.add_subplot(2, 2, 1)
    _plot_umap(adata, ax=sub, color='S_score')
    # sub.set_aspect('equal')
    hf.set_axes_label(x_label='UMAP_1', y_label='UMAP_2', axes=sub)
    sub.set_title('S-Phase cluster')

    sub2 = fig_2d_visu.add_subplot(2, 2, 2)
    _plot_umap(adata, ax=sub2, color='G2M_score')
    # sub2.set_aspect('equal')
    sub2.set_title('G2- and M-Phase clusters')
    hf.set_axes_label(x_label='UMAP_1', y_label='UMAP_2', axes=sub2)

    sub3 = fig_2d_visu.add_subplot(2, 2, 3)
    _plot_umap(adata, ax=sub3, color='phase')
    # sub3.set_aspect('equal')
    sub3.set_title('Phase clusters')
    hf.set_axes_label(x_label='UMAP_1', y_label='UMAP_2', axes=sub3)

    sub4 = fig_2d_visu.add_subplot(2, 2, 4)
    _plot_umap(adata, ax=sub4, color='n_counts')
    # sub4.set_aspect('equal')
    sub4.set_title('2D counts')
    hf.set_axes_label(x_label='UMAP_1', y_label='UMAP_2', axes=sub4)

    hf.set_title(sup_title="2D Cell Cycle Visualization", figure=fig_2d_visu)
    fig_2d_visu.tight_layout(pad=1.5, w_pad=1.5, h_pad=1.5)
    fig_2d_visu.subplots_adjust(top=0.88)
    # fig_2d_visu.gca().set_aspect('equal', adjustable='box')

    if raw:
        plt.savefig(os.path.join(save_folder, "raw_2D_Cell_Cycle.png"))
    else:
        plt.savefig(os.path.join(save_folder, "2D_Cell_Cycle.png"))
    plt.close()


def plot_highest_expr_genes(adata, save_folder, raw=False):
    """Plot top 20 highly expressed genes

    Parameters
    ----------
    adata : annData
    save_folder : str
    raw : bool

    Returns
    -------

    """

    # whole data set
    fig_2d_visu = plt.figure(facecolor='w', edgecolor='k', figsize=fig_size)
    sub = fig_2d_visu.add_subplot(1, 1, 1)
    sc.pl.highest_expr_genes(adata, n_top=20, ax=sub, show=False)
    hf.set_title(sup_title="Highest expressed Genes", figure=fig_2d_visu)
    fig_2d_visu.subplots_adjust(bottom=0.125, left=0.125, top=0.875, right=0.875)

    if raw:
        plt.savefig(os.path.join(save_folder, "whole_data_set_raw_HEG.png"))
    else:
        plt.savefig(os.path.join(save_folder, "whole_data_set_HEG.png"))

    # per sample
    sample_names = np.unique(adata.obs['sample'])
    for sample in sample_names:
        adata_sample = adata.copy()
        adata_sample = adata_sample[adata_sample.obs["sample"] == sample]
        fig_2d_visu = plt.figure(facecolor='w', edgecolor='k', figsize=fig_size)
        sub = fig_2d_visu.add_subplot(1, 1, 1)
        sc.pl.highest_expr_genes(adata_sample, n_top=20, ax=sub, show=False)
        hf.set_title(sup_title="Highest expressed Genes", figure=fig_2d_visu)
        fig_2d_visu.subplots_adjust(bottom=0.125, left=0.125, top=0.875, right=0.875)

        if raw:
            plt.savefig(os.path.join(save_folder, "_".join([sample, "raw_HEG.png"])))
        else:
            plt.savefig(os.path.join(save_folder, "_".join([sample, "HEG.png"])))
    plt.close()


def plot_highly_variable_genes(adata, type_dataset, save_folder, raw=False):
    """Plot highly variable genes (HVG)

    Parameters
    ----------
    adata : annData
    type_dataset : str
    save_folder : str
    raw : bool

    Returns
    -------

    """
    # whole data set
    sc.pl.highly_variable_genes(adata, show=False)

    if raw:
        plt.savefig(os.path.join(save_folder, "_".join([type_dataset, "raw_HVG.png"])))
    else:
        plt.savefig(os.path.join(save_folder, "_".join([type_dataset, "HVG.png"])))

    plt.close()


def plot_pc_combs(adata, type_dataset, save_folder, raw=False):
    """Plot explained variance by Principal Components (PCs)
    Determine the number of informative principal components
    ==> look at variance contribution of each component
    ==> apply Elbow method

    Parameters
    ----------
    adata : annData
    type_dataset : str
        e.g. if batch corrected or not
    save_folder : str
    raw : bool
        if raw annData is provided

    Returns
    -------

    """
    sc.pl.pca_variance_ratio(adata, show=False)

    if raw:
        plt.savefig(os.path.join(save_folder, "{}_raw_pc_combs.png".format(type_dataset)))
    else:
        plt.savefig(os.path.join(save_folder, "{}_pc_combs.png".format(type_dataset)))
    plt.close()


def plot_batch_correction(adata, save_folder, batch_key, possible_batch_effect):
    """Plot batch correction results

    Parameters
    ----------
    adata : annData
    save_folder : str
    batch_key : str
        kind of data matrix used for visualisation
    possible_batch_effect : str

    Returns
    -------

    """

    plot_scatter_contour(adata=adata, batch_key=batch_key, observable=possible_batch_effect, save_folder=save_folder,
                         show_points=True)
    fig = plt.figure(facecolor='w', edgecolor='k', figsize=(14, 8))
    sub_1 = fig.add_subplot(1, 1, 1)
    sc.pl.umap(adata, color=possible_batch_effect, ax=sub_1, show=False, title="Sample")
    fig.subplots_adjust(bottom=0.125, left=0.025, top=0.875, right=0.850)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "scanorama_{}_{}.png".format(batch_key, possible_batch_effect)))
    plt.close()


def plot_scatter_contour(adata, batch_key, observable, save_folder, show_points=True):
    """Draw contours around clusters
    Author: Christina K. L. Hillig

    Parameters
    ----------
    adata : annData
    batch_key : str
    observable : str
        column in adata.obs
    save_folder : str
    show_points : bool

    Returns
    -------

    """
    coordinates = []
    for spot_counter in range(adata.shape[0]):
        row = [adata.obs[observable][spot_counter],
               adata.obsm['X_umap'][spot_counter][0], adata.obsm['X_umap'][spot_counter][1]]
        coordinates.append(row)

    df = pd.DataFrame(coordinates, columns=[observable, 'UMAP 1', 'UMAP 2'])
    df = df.sort_values(by=[observable], ascending=True)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot()
    sns.set_style("white")

    x = 'UMAP 1'
    y = 'UMAP 2'
    variables = np.unique(df[observable].values)
    print(sorted(variables))

    if show_points:
        if 10 < len(variables) <= 20:
            for label, color in zip(sorted(np.unique(df[observable].values)),
                                    sns.color_palette("tab20")[:len(variables)]):
                temp_df = df.loc[df[observable] == label]
                sns.kdeplot(temp_df[x], temp_df[y], color=color, shade_lowest=False, alpha=0.1, ax=ax)
                sns.kdeplot(temp_df[x], temp_df[y], shade=True, shade_lowest=False, alpha=0.5, color=color, ax=ax)
            sns.scatterplot(x=x, y=y, hue=observable, data=df, palette=sns.color_palette("tab20")[:len(variables)],
                            ax=ax, s=5)
        elif 20 < len(variables) < 100:
            for label, color in zip(sorted(np.unique(df[observable].values)),
                                    sc.pl.palettes.godsnot_102[:len(variables)]):
                temp_df = df.loc[df[observable] == label]
                sns.kdeplot(temp_df[x], temp_df[y], color=color, shade_lowest=False, alpha=0.1, ax=ax)
                sns.kdeplot(temp_df[x], temp_df[y], shade=True, shade_lowest=False, alpha=0.5, color=color, ax=ax)
            sns.scatterplot(x=x, y=y, hue=observable, data=df, palette=sc.pl.palettes.godsnot_102[:len(variables)],
                            ax=ax, s=5)
        else:
            for label, color in zip(sorted(np.unique(df[observable].values)),
                                    sns.color_palette("tab10")[:len(variables)]):
                temp_df = df.loc[df[observable] == label]
                sns.kdeplot(temp_df[x], temp_df[y], color=color, shade_lowest=False, alpha=0.1, ax=ax)
                sns.kdeplot(temp_df[x], temp_df[y], shade=True, shade_lowest=False, alpha=0.5, color=color, ax=ax)
            sns.scatterplot(x=x, y=y, hue=observable, data=df, palette=sns.color_palette("tab10")[:len(variables)],
                            ax=ax, s=5)
    else:
        if 10 < len(variables) <= 20:
            for label, color in zip(sorted(np.unique(df[observable].values)),
                                    sns.color_palette("tab20")[:len(variables)]):
                temp_df = df.loc[df[observable] == label]
                sns.kdeplot(temp_df[x], temp_df[y], color=color, shade_lowest=False, alpha=0.1, ax=ax)
                sns.kdeplot(temp_df[x], temp_df[y], shade=True, shade_lowest=False, alpha=0.5, color=color, ax=ax)
        elif 20 < len(variables) < 100:
            for label, color in zip(sorted(np.unique(df[observable].values)),
                                    sc.pl.palettes.godsnot_102[:len(variables)]):
                temp_df = df.loc[df[observable] == label]
                sns.kdeplot(temp_df[x], temp_df[y], color=color, shade_lowest=False, alpha=0.1, ax=ax)
                sns.kdeplot(temp_df[x], temp_df[y], shade=True, shade_lowest=False, alpha=0.5, color=color, ax=ax)
        else:
            for label, color in zip(sorted(np.unique(df[observable].values)),
                                    sns.color_palette("tab10")[:len(variables)]):
                temp_df = df.loc[df[observable] == label]
                sns.kdeplot(temp_df[x], temp_df[y], color=color, shade_lowest=False, alpha=0.1, ax=ax)
                sns.kdeplot(temp_df[x], temp_df[y], shade=True, shade_lowest=False, alpha=0.5, color=color, ax=ax)

    # Put the legend out of the figure
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()

    if save_folder is None:
        plt.show()
    else:
        plt.savefig(os.path.join(save_folder, "_".join([observable, "UMAP_{}_counterplot.png".format(batch_key)])))
    plt.close()
