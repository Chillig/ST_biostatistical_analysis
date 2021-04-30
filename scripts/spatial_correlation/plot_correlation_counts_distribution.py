"""Visualise Weighted Spatial Pearson Correlation and Connectivity graphs with and without tissue image
    File name: plot_correlation_counts_distribution.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

# import scripts
from scripts.spatial_correlation import corr_statistics as corr_stats
from scripts.spatial_correlation import helper_functions

# Plotting packages
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import collections as mc
import seaborn as sns

# System specific
import os

# Calculation packages
from collections import OrderedDict
import scanpy as sc
import numpy as np
import pandas as pd
import pingouin
import scipy.stats as scstats

# Figure params
sc.set_figure_params(color_map='viridis')
fig_size, title_fontsize, axis_label_fontsize, legend_fontsize, fileformat, img_key, xy_ticks, text_fontsize = \
    helper_functions.figure_params()


def plot_clusters_counts(sup_adata, cyto, save_folder, obs_label, obs_counts):
    """

    :param sup_adata: [dict] containing annData object
    :param cyto: [string] cytokine name
    :param save_folder: [string]
    :param obs_label: [string]
    :param obs_counts: [string]
    :return:
    """
    # Size of count spots
    size = 0.9

    color_dict = OrderedDict()
    if cyto == 'IFNG':
        color_dict[cyto] = "#ff7f00"  # orange LICHEN "#ff7f00"
    elif cyto == "IL13":
        color_dict[cyto] = "#e41a1c"  # red AE
    else:
        color_dict[cyto] = "#377eb8"  # blue PSO

    color_dict["Responders"] = "y"
    # color_dict[" & ".join([cyto, "Responders"])] = "gold"
    color_dict["Others"] = "silver"

    diseases = ["PSO", "AE", "LICHEN", "PRP"]
    biopsy_type = ["LESONAL", "NON LESIONAL"]

    # samples = np.unique(sup_adata[cyto].obs['sample'])
    samples, crops_img = helper_functions.get_cropped_sampleimg(img_key=img_key)
    for ind, sample in enumerate(samples):
        temp_adata = sup_adata[sup_adata.obs["sample"] == sample]
        # read out counts of cytokine and responder spots
        test_counts = temp_adata[temp_adata.obs[obs_label] != "Others"]
        # get labels of cytokine and responder spots
        cell_types_unique = list(np.unique(temp_adata.obs[obs_label]))

        list_colors = []
        for i_color in cell_types_unique:
            list_colors.append(color_dict[i_color])

        if len(cell_types_unique) > 0:

            diagnose = temp_adata.obs[diseases].loc[:, (temp_adata.obs[diseases] != 0).all()].columns.values[0]
            fig, ax = plt.subplots(figsize=fig_size)
            sc.pl.spatial(temp_adata, size=1.3, img_key=img_key, library_id=sample,
                          color=obs_label, groups=cell_types_unique, zorder=1, show=False,
                          legend_loc='left margin', ax=ax, alpha_img=1, crop_coord=crops_img[ind],
                          palette=list_colors)
            if cyto == "IL17A":
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=sample,
                              color=obs_counts, zorder=2, alpha_img=0, show=False, crop_coord=crops_img[ind],
                              ax=ax, vmin=0, vmax=800,
                              title=" ".join(["Diagnose:", diagnose, "; Biopsy type:",
                                              temp_adata.obs[biopsy_type].loc[:, (temp_adata.obs[biopsy_type] !=
                                                                                  0).all()].columns.values[0]]))
            elif cyto == "IFNG":
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=sample,
                              color=obs_counts, zorder=2, alpha_img=0, show=False, crop_coord=crops_img[ind],
                              ax=ax, vmin=0, vmax=100,
                              title=" ".join(["Diagnose:", diagnose, "; Biopsy type:",
                                              temp_adata.obs[biopsy_type].loc[:, (temp_adata.obs[biopsy_type] !=
                                                                                  0).all()].columns.values[0]]))
            elif cyto == "IL13":
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=sample,
                              color=obs_counts, zorder=2, alpha_img=0, show=False, crop_coord=crops_img[ind],
                              ax=ax, vmin=0, vmax=50,
                              title=" ".join(["Diagnose:", diagnose, "; Biopsy type:",
                                              temp_adata.obs[biopsy_type].loc[:, (temp_adata.obs[biopsy_type] !=
                                                                                  0).all()].columns.values[0]]))
            else:
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=sample,
                              color=obs_counts, zorder=2, alpha_img=0, show=False, crop_coord=crops_img[ind],
                              ax=ax, vmin=0, vmax=100,
                              title=" ".join(["Diagnose:", diagnose, "; Biopsy type:",
                                              temp_adata.obs[biopsy_type].loc[:, (temp_adata.obs[biopsy_type] !=
                                                                                  0).all()].columns.values[0]]))
            # Invert both axis due to flipped and mirrored images
            ax.invert_xaxis()
            ax.invert_yaxis()

            # Colorbar label
            cbar_ax = fig.axes[-1]
            cbar_ax.get_yaxis().labelpad = 15
            cbar_ax.set_ylabel('UMI-counts', rotation=90, fontsize=legend_fontsize)
            cbar_ax.tick_params(labelsize=xy_ticks)

            # Legend: outside of axis
            leg = ax.legend(bbox_to_anchor=(2, 0.5), ncol=1, fontsize=legend_fontsize)
            leg.get_frame().set_linewidth(0.0)

            plt.tight_layout()

            if img_key == "lowres":
                plt.savefig(os.path.join(save_folder, "_".join(["Radial_plot", sample, cyto, ".png"])),
                            bbox_inches='tight',  bbox_extra_artists=(leg,))
                plt.close()
            else:
                plt.savefig(os.path.join(save_folder, "_".join(["Radial_plot", sample, cyto, fileformat])),
                            bbox_inches='tight',  bbox_extra_artists=(leg,))
                plt.close()


def plot_nn_spots(adata, sample, obs, save_folder, title, distance=1):
    """
    Plot nearest neighbor spot and center spot

    :param adata: [annData]
    :param sample: [string]
    :param obs: [string] observation group to plot
    :param save_folder: [string]
    :param title: [string]
    :param distance: [int]
    :return:
    """
    # Create new folder in this directory for specific distance
    save_folder = os.path.join(save_folder, str(distance))
    os.makedirs(save_folder, exist_ok=True)

    # Size of count spots
    size = 0.9
    fig, ax = plt.subplots(figsize=fig_size)
    sc.pl.spatial(adata, size=size, img_key='hires', library_id=sample, color=obs, title=title, show=False, ax=ax,
                  color_map='viridis')
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "_".join(["ST_NN_spots", title, fileformat])))
    plt.close()


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


def plot_normedcounts(dfx, dfy, index_counter, index_counter_sample, cyto, sample, save_folder, distance=1):
    """Plot normalised counts of cytokines and responders

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
    ax.set_xlabel('Normed responder genes', fontsize=axis_label_fontsize)
    ax.set_ylabel(" ".join(["Normed", cyto, 'counts']), fontsize=axis_label_fontsize)
    ax.set_title(' '.join(['Distance', str(distance), sample, cyto]), fontsize=title_fontsize)
    ax.set_yticks(np.arange(0, np.amax(dfy[index_counter - index_counter_sample:index_counter]) + 1, 1))
    ax.tick_params(labelsize=xy_ticks)
    plt.tight_layout()
    fig.savefig(
        os.path.join(save_folder, '_'.join(['Distance', str(distance), sample, cyto, "Normed", fileformat])))
    plt.close()


def plot_st_correlation(df_counts, cytokine_responders, save_folder, distance):
    """Calculate spatial Correlation between each cytokine positive spot and its nn responder genes spots for
        each cluster

    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    distance : int

    Returns
    -------

    """

    pvals = []
    for cyto in cytokine_responders:
        resp_name = "_".join([cyto, 'responder'])
        temp_df = df_counts[[cyto, resp_name, "Cluster_size", "disease"]].copy()

        # temp_df = temp_df.replace({np.nan: 0})
        temp_df = temp_df[~np.isnan(temp_df[cyto].values.astype(np.float64))]
        temp_df[cyto] = temp_df[cyto].values.astype(np.float)
        temp_df["_".join([cyto, 'responder'])] = temp_df["_".join([cyto, 'responder'])].values.astype(np.float)

        # 1. Calculate correlation
        # check if data is normally distributed
        p_w = corr_stats.apply_wilcoxontest(df_highcounts=temp_df["_".join([cyto, 'responder'])],
                                            df_lowcounts=temp_df[cyto])
        # if p_w:
        # Always report the correlation and the p-value:
        # -> Use pearson correlation and spearman and compare them
        # -> at low statisitcs the p-value might be infiltrated
        sig_r = pingouin.corr(x=temp_df["_".join([cyto, 'responder'])],
                              y=temp_df[cyto], method='pearson')
        pvals.append(sig_r['p-val'].values[0])

        get_sizes = np.unique(temp_df['Cluster_size'])
        labels_names = ['fit', 'confidence interval']
        labels_names.extend(get_sizes)

        # Increase point size
        temp_df["Cluster size"] = (temp_df["Cluster_size"].values**1).astype(int)

        #  Plot correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        sns.regplot(data=temp_df, x=resp_name, y=cyto, ax=ax, scatter=False, color="black", label=None)
        if len(np.unique(temp_df['Cluster_size'])) > 1:
            sns.scatterplot(data=temp_df, x=resp_name, y=cyto, size="Cluster size", legend='full',
                            hue="Cluster size", ax=ax,
                            palette=sns.color_palette("husl")[:len(np.unique(temp_df['Cluster_size']))])
        else:
            sns.scatterplot(data=temp_df, x=resp_name, y=cyto, size="Cluster size", legend=False,
                            hue="Cluster size", ax=ax, palette=['k'])

        # if p_w:
        ax.text(temp_df.max()[1] / 2 - temp_df.max()[1]/10, temp_df.max()[0] + temp_df.max()[0] / 8,
                'r = {:.2f}; p = {:.2e}'.format(sig_r['r'].values[0], sig_r['p-val'].values[0]),
                fontstyle='italic', fontsize=text_fontsize)
        ax.set_xlabel(" ".join(["Responder Counts"]))
        if cyto == 'IFNG':
            ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel(" ".join([cyto, 'Counts']), fontsize=axis_label_fontsize)
        ax.set_yticks(np.arange(0, temp_df.max()[0] + 1, 1))
        ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 50])
        ax.set_ylim([-0.5, temp_df.max()[0] + temp_df.max()[0] / 5])
        ax.tick_params(labelsize=xy_ticks)

        plt.tight_layout()
        sns.despine(ax=ax)
        fig.savefig(os.path.join(save_folder, "_".join(['Fig4A', str(distance), cyto, resp_name, fileformat])))
        plt.close()

    return pvals


def plot_spatial_correlation(df_counts, cytokine_responders, save_folder, distance):
    """Calculate spatial Correlation between each Cyto+ spot and its nn responder genes spots for each cluster

    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    distance : int

    Returns
    -------
    list of p-values

    """

    p_vals = []
    for cyto in cytokine_responders:
        resp_name = "_".join([cyto, 'responder'])
        temp_df = df_counts[[cyto, resp_name]].copy()

        temp_df = temp_df[~np.isnan(temp_df[cyto].values.astype(np.float64))]
        temp_df[cyto] = temp_df[cyto].values.astype(np.float)
        temp_df["_".join([cyto, 'responder'])] = temp_df["_".join([cyto, 'responder'])].values.astype(np.float)

        # 1. Calculate correlation
        # 1.1 Apply weights to cytokine counts and add it to df: Multiply with cluster size
        # Always report the correlation and the p-value:
        # -> Use pearson correlation and spearman and compare them
        # -> at low statisitcs the p-value might be infiltrated
        sig_r = pingouin.corr(x=temp_df["_".join([cyto, 'responder'])],
                              y=temp_df[cyto], method='pearson')
        p_vals.append(sig_r['p-val'].values[0])

        #  2. Plot correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        sns.regplot(data=temp_df, x=resp_name, y=cyto, ax=ax, scatter=False, color="black", label=None)
        ax.scatter(data=temp_df, x=resp_name, y=cyto, c='k')
        # Add text: Correlation value and p-value
        # Always report the correlation and the p-value:
        # -> Use pearson correlation and spearman and compare them
        # -> at low statistics the p-value might be infiltrated
        ax.text(temp_df.max()[1] / 2 - temp_df.max()[1]/10, temp_df.max()[0],
                'r = {:.2f}; p = {:.2e}'.format(sig_r['r'].values[0], sig_r['p-val'].values[0]),
                fontstyle='italic', fontsize=text_fontsize)

        # Axis params
        ax.set_xlabel(" ".join(["Responder Counts"]), fontsize=axis_label_fontsize)
        if cyto == 'IFNG':
            ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel(" ".join([cyto, 'Counts']), fontsize=axis_label_fontsize)
        # ax.set_yticks(np.arange(0, temp_df.max()[0] + 1, 1))
        ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 50])
        ax.set_ylim([-0.5, temp_df.max()[0] + temp_df.max()[0] / 50])
        ax.tick_params(labelsize=xy_ticks)

        plt.tight_layout()
        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure
        fig.savefig(os.path.join(save_folder, "_".join(['Fig1G', str(distance), cyto, resp_name, fileformat])))
        plt.close()

    return p_vals


def plot_st_weightedcorrelation(df_counts, cytokine_responders, save_folder, distance):
    """Calculate spatial, weighted Correlation between each Cyto+ spot and its nn responder genes spots for each cluster

    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    distance : int

    Returns
    -------
    list of p-values

    """
    # Get unique combiations of the tissue layer to assign to each a color before the plotting
    # 1. map each string to an integer value
    color_seq = pd.factorize(df_counts['tissue_layer'])[0]
    num_unique_colors = len(np.unique(color_seq))
    # 2. assign a color to each combination
    dict_tissue_color = dict(zip(df_counts['tissue_layer'], color_seq))
    tissuecomb_colors = plt.cm.get_cmap("tab20", num_unique_colors)

    p_vals = []
    for cyto in cytokine_responders:
        resp_name = "_".join([cyto, 'responder'])
        temp_df = df_counts[[cyto, resp_name, "Cluster_size", "disease", 'tissue_layer']].copy()
        weighted_cytoname = "_".join(['weighted', cyto])

        temp_df = temp_df[~np.isnan(temp_df[cyto].values.astype(np.float64))]
        temp_df[cyto] = temp_df[cyto].values.astype(np.float)
        temp_df["_".join([cyto, 'responder'])] = temp_df["_".join([cyto, 'responder'])].values.astype(np.float)

        # # map each string to an integer value
        cyto_colorseq = []
        for ind_dict, value in enumerate(temp_df['tissue_layer']):
            cyto_colorseq.append(dict_tissue_color[value])
        num_colors = len(np.unique(cyto_colorseq))

        # 1. Calculate correlation
        # 1.1 Apply weights to cytokine counts and add it to df: Multiply with cluster size
        weighted_cytocounts = temp_df[cyto] * temp_df['Cluster_size']
        temp_df[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)
        # 1.2 Calculate weighted Pearson Correlation
        correlation = corr_stats.weighted_corr(x=temp_df["_".join([cyto, 'responder'])], y=temp_df[cyto],
                                               w=temp_df['Cluster_size'])  # np.ones(len(temp_df.index))
        # 1.3 Calculate p-value
        p_w = corr_stats.apply_statstest(df_highcounts=temp_df["_".join([cyto, 'responder'])],
                                         df_lowcounts=weighted_cytocounts, correlation_value=correlation)

        p_vals.append(p_w)

        get_sizes = np.unique(temp_df['Cluster_size'])

        # Increase point size
        temp_df["Cluster size"] = (temp_df["Cluster_size"].values**3 + 10).astype(int)

        #  2. Plot correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        sns.regplot(data=temp_df, x=resp_name, y=cyto, ax=ax, scatter=False, color="black", label=None)
        sc_pl = ax.scatter(data=temp_df, x=resp_name, y=cyto, s="Cluster size", c=cyto_colorseq,
                           cmap=ListedColormap(tissuecomb_colors.colors[np.unique(cyto_colorseq)]))
        # Add text: Correlation value and p-value
        # if p_w:
        # Always report the correlation and the p-value:
        # -> Use pearson correlation and spearman and compare them
        # -> at low statistics the p-value might be infiltrated
        ax.text(temp_df.max()[1] / 2 - temp_df.max()[1]/10, temp_df.max()[0] + temp_df.max()[0] / 8,
                'r = {:.2f}; p = {:.2e}'.format(correlation, p_w), fontstyle='italic', fontsize=text_fontsize)

        # Axis params
        ax.set_xlabel(" ".join(["Responder Counts"]), fontsize=axis_label_fontsize)
        if cyto == 'IFNG':
            ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel(" ".join([cyto, 'Counts']), fontsize=axis_label_fontsize)
        ax.set_yticks(np.arange(0, temp_df.max()[0] + 1, 1))
        ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 50])
        ax.set_ylim([-0.5, temp_df.max()[0] + temp_df.max()[0] / 5])
        ax.tick_params(labelsize=xy_ticks)

        # Make a legend for cytokine cluster size
        for ind, pw in enumerate(np.unique(temp_df["Cluster size"].values)):
            if pw != 0:
                ax.scatter([], [], s=pw, c="k", label=" {}".format(get_sizes[ind]))

        legend_handles, legend_labels = plt.gca().get_legend_handles_labels()
        if any(get_sizes > 9):
            label_space = 2
        else:
            label_space = 0.8
        plt.legend(legend_handles[1:], legend_labels[1:], labelspacing=label_space, title="Cluster size",
                   title_fontsize=title_fontsize,
                   fontsize=legend_fontsize,  # size of title
                   bbox_transform=fig.transFigure, loc='best',  # bbox_to_anchor=(1.01, 1),  # place legend outside
                   borderpad=1, scatterpoints=1, handletextpad=2, handlelength=2,
                   frameon=True, framealpha=0.6, facecolor="w")

        # Colorbar to add tissue layers as colors
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('bottom', size='10%', pad=0.6)
        cbar = fig.colorbar(sc_pl, cax=cax, orientation='horizontal')
        cax.set_xlabel('Tissue layers', fontsize=axis_label_fontsize)
        # cbar.set_ticks(np.linspace(0, num_colors - 1, num_colors + 1)[:-1] + 0.5)
        cbar.set_ticks(np.linspace(0, num_colors - (num_colors - np.array(cyto_colorseq).max()), num_colors + 1)[:-1]
                       + num_colors / (np.array(cyto_colorseq).max() - np.array(cyto_colorseq).min()))
        cbar.ax.set_xticklabels(list(temp_df['tissue_layer'].unique()), fontsize=legend_fontsize, rotation=-45,
                                ha="left")

        plt.tight_layout()
        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure
        fig.savefig(os.path.join(save_folder, "_".join(['Fig4C', str(distance), cyto, resp_name, fileformat])))
        plt.close()

    return p_vals


def plot_weighted_spatialcorrelation_wotissuelabels(df_counts, cytokine_responders, save_folder, distance):
    """Calculate Correlation between each Cytokine positive spot and its nn responder genes spots for each sample

    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    distance : int

    Returns
    -------

    """
    # Get unique combiations of the tissue layer to assign to each a color before the plotting
    # 1. map each string to an integer value
    color_seq = pd.factorize(df_counts['tissue_layer'])[0]
    num_unique_colors = len(np.unique(color_seq))
    # 2. assign a color to each combination
    dict_tissue_color = dict(zip(df_counts['tissue_layer'], color_seq))
    tissuecomb_colors = plt.cm.get_cmap("tab20", num_unique_colors)

    # 3. Plot Colorbar with tissue layer combinations
    plot_standalone_colorbar(tissuecomb_colors=tissuecomb_colors, labels=np.unique(df_counts['tissue_layer']),
                             save_folder=save_folder)

    # 4. get unique cluster sizes
    cluster_sizes = np.unique(df_counts['Cluster_size'])
    plot_standalone_legend(cluster_sizes, save_folder=save_folder)

    p_vals = []
    for cyto in cytokine_responders:
        resp_name = "_".join([cyto, 'responder'])
        temp_df = df_counts[[cyto, resp_name, "Cluster_size", "disease", 'tissue_layer']].copy()
        weighted_cytoname = "_".join(['weighted', cyto])

        temp_df = temp_df[~np.isnan(temp_df[cyto].values.astype(np.float64))]
        temp_df[cyto] = temp_df[cyto].values.astype(np.float)
        temp_df["_".join([cyto, 'responder'])] = temp_df["_".join([cyto, 'responder'])].values.astype(np.float)

        # # map each string to an integer value
        cyto_colorseq = []
        for ind_dict, value in enumerate(temp_df['tissue_layer']):
            cyto_colorseq.append(dict_tissue_color[value])

        # 1. Calculate correlation
        # 1.1 Apply weights to cytokine counts and add it to df: Multiply with cluster size
        weighted_cytocounts = temp_df[cyto] * temp_df['Cluster_size']
        temp_df[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)
        # 1.2 Calculate weighted Pearson Correlation
        correlation = corr_stats.weighted_corr(
            x=temp_df["_".join([cyto, 'responder'])], y=temp_df[cyto], w=temp_df['Cluster_size'])
        # 1.3 Calculate p-value
        p_w = corr_stats.apply_statstest(df_highcounts=temp_df["_".join([cyto, 'responder'])],
                                         df_lowcounts=weighted_cytocounts, correlation_value=correlation)
        p_vals.append(p_w)

        get_sizes = np.unique(temp_df['Cluster_size'])

        # Increase point size
        temp_df["Cluster size"] = (temp_df["Cluster_size"].values**2 * 36).astype(int)

        #  2. Plot correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        ax_sns = sns.regplot(data=temp_df, x=resp_name, y=cyto, ax=ax, scatter=False, color="black", label=None)
        plt.setp(ax_sns.lines, zorder=1)
        sc_pl = ax.scatter(data=temp_df, x=resp_name, y=cyto, s="Cluster size", c=cyto_colorseq,
                           cmap=ListedColormap(tissuecomb_colors.colors[np.unique(cyto_colorseq)]), zorder=2)

        # Axis params
        ax.set_xlabel(" ".join(["Responder Counts"]), fontsize=axis_label_fontsize)
        if cyto == 'IFNG':
            ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel(" ".join([cyto, 'Counts']), fontsize=axis_label_fontsize)
        if any(get_sizes > 6):
            # ax.set_yticks(np.arange(0, temp_df.max()[0] + 3, 1))
            # if temp_df.max()[0] < 10:
            #     ax.set_yticks(np.arange(0, temp_df.max()[0] + 1, 1))
            # else:
            #     ax.set_yticks(np.arange(0, temp_df.max()[0] + 3, 2))
            # ax.set_ylim([-0.5, temp_df.max()[0] + 2])
            # ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 13])
            # Add text: Correlation value and p-value
            ax.text(temp_df.max()[1] / 2 - temp_df.max()[1] / 10, temp_df.max()[0] + 1,
                    'r = {:.2f}; p = {:.2e}'.format(correlation, p_w), fontstyle='italic', fontsize=text_fontsize)

        else:
            # ax.set_yticks(np.arange(0, temp_df.max()[0] + 1, 1))
            # ax.set_ylim([-0.5, temp_df.max()[0] + temp_df.max()[0] / 20])
            # ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 20])
            # Add text: Correlation value and p-value
            ax.text(temp_df.max()[1] / 2 - temp_df.max()[1] / 10, temp_df.max()[0],
                    'r = {:.2f}; p = {:.2e}'.format(correlation, p_w), fontstyle='italic', fontsize=text_fontsize)
        ax.tick_params(labelsize=xy_ticks)

        # # Make a legend for cytokine cluster size
        # for ind, pw in enumerate(np.unique(temp_df["Cluster size"].values)):
        #     if pw != 0:
        #         ax.scatter([], [], s=pw, c="k", label=" {}".format(get_sizes[ind]))
        # legend_handles, legend_labels = plt.gca().get_legend_handles_labels()
        # if any(get_sizes > 9):
        #     label_space = 2
        # else:
        #     label_space = 0.8

        # leg = ax.legend(legend_handles[1:], legend_labels[1:], labelspacing=label_space, title="Cluster size",
        #                  title_fontsize=title_fontsize, fontsize=legend_fontsize,  # size of title
        #                  bbox_transform=fig.transFigure, loc='best',
        #                  borderpad=1, scatterpoints=1, handletextpad=2, handlelength=2,
        #                  bbox_to_anchor=(1., 1.),
        #                  frameon=True, framealpha=0.6, facecolor="w")

        # # Colorbar to add tissue layers as colors but dont show the labels in the colorbarticks
        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes('bottom', size='10%', pad=0.6)
        # cbar = fig.colorbar(sc_pl, cax=cax, orientation='horizontal')
        # # Remove ticks and labels
        # cbar.set_ticks([])
        # cax.set_xlabel('Tissue layers', fontsize=axis_label_fontsize)

        plt.tight_layout()

        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure
        fig.savefig(os.path.join(save_folder, "_".join(['Fig4A', str(distance), cyto, resp_name, fileformat])),)
                    # bbox_inches='tight',  bbox_extra_artists=(leg,))
        plt.close()

    return p_vals


def plot_standalone_colorbar(tissuecomb_colors, labels, save_folder):
    """Plot a standalone vertical and horizontal Colorbar

    Parameters
    ----------
    tissuecomb_colors : matplotlib.colors.ListedColormap
    labels : list of str
    save_folder : str

    Returns
    -------

    """

    num_colors = len(tissuecomb_colors.colors)
    norm = mpl.colors.Normalize(vmin=0, vmax=num_colors)
    # replace _ with  &
    labels = [label.replace("_", " & ") for label in labels]

    fig = plt.figure(figsize=fig_size)
    # [left, bottom, width, height]
    ax = fig.add_axes([0.1, 0.05, 0.16, 0.9])  # work only if rotation = 0
    cbar = mpl.colorbar.ColorbarBase(ax, cmap=ListedColormap(tissuecomb_colors.colors), norm=norm,
                                     orientation='vertical')
    cbar.set_ticks(np.linspace(0, num_colors, num_colors + 1)[:-1] + 0.5)
    cbar.ax.set_yticklabels(labels, fontsize=xy_ticks, rotation=0, ha="left")
    ax.set_title('Tissue layers', fontsize=title_fontsize)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "colorbar.pdf"))
    plt.close()


def plot_standalone_legend(points, save_folder):
    """Plot a legend without a graph

    Parameters
    ----------
    points : numpy.array
    save_folder : str

    Returns
    -------

    """
    # increase point size
    large_points = points.astype(int)**2 * 36
    colors = ['k'] * len(points)

    # Draw circles
    patches = [plt.scatter([], [], marker="o", s=large_points.astype(int)[i], color=colors[i],
                           label="{:s}".format(str(points[i]))) for i in range(len(points))]
    plt.close()

    if any(points > 6):
        label_space = 3
    else:
        label_space = 1

    # Create figure without an axis
    fig = plt.figure(figsize=fig_size)
    fig.legend(patches, points, labelspacing=label_space, title="Cluster size",
               loc='center', frameon=False, facecolor="w", handletextpad=2, handlelength=2,
               title_fontsize=title_fontsize, fontsize=legend_fontsize,
               bbox_transform=fig.transFigure, borderpad=.0, scatterpoints=1)
    # Save figure
    fig.savefig(os.path.join(save_folder, "_".join(["Legend", fileformat])),)
    plt.close()


def plot_evaluate_distance(significance, cytokines, save_folder):
    """Elbow plot for best distance/ radius evaluation

    Parameters
    ----------
    significance : list
    cytokines : list
    save_folder : str

    Returns
    -------

    """
    # Evaluate distance via elbow plot
    significance = np.array(significance).T
    for ind, cyto in enumerate(cytokines):
        fig, ax = plt.subplots(figsize=fig_size)
        ax.plot(np.arange(1, len(significance.T) + 1), -np.log10(significance[ind, :]))
        ax.set_xlabel('Radius', fontsize=axis_label_fontsize)
        ax.set_ylabel(r'-log$_{10}$(p-values)', fontsize=axis_label_fontsize)
        ax.set_xticks(np.arange(1, len(significance.T) + 1))
        ax.set_title(cyto)

        plt.tight_layout()
        fig.savefig(os.path.join(save_folder, "_".join(['Radius_Evaluation', cyto, fileformat])))
        plt.close()
        print("Highest significance at radius: ", np.argmax(-np.log10(significance[ind, :])) + 1)


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


def plot_graphs(all_coordinates, cluster_coordinates, graph, save_folder, title):
    """Plot points in coordiante system and mark cyto+ spot in red and responder spots in blue
        Connect nodes of cyto+ spots

    Parameters
    ----------
    all_coordinates : list of tuples
    cluster_coordinates : list of tuples
    graph : list of graphs
    save_folder : str
    title : str

    Returns
    -------

    """
    # Get edges to draw between nodes (cyto+ spots)
    edges = []
    for ind_i, i in enumerate(graph):
        edges.append([])
        for j in i:
            edges[ind_i].append([cluster_coordinates[j][0], cluster_coordinates[j][1]])

    # edges: [[[x1, y1], [x2, y2]], [[x3, y3], [x4, y4]], ..]

    max_xlim = np.amax((np.array(all_coordinates)[:, 0]))
    max_ylim = np.amax((np.array(all_coordinates)[:, 1]))

    min_xlim = np.amin((np.array(all_coordinates)[:, 0]))
    min_ylim = np.amin((np.array(all_coordinates)[:, 1]))

    delta_x = (max_xlim - min_xlim)
    delta_y = (max_ylim - min_ylim)

    greates_interval = np.amax([delta_x, delta_y])

    # Plotting example
    fig, ax = plt.subplots(figsize=fig_size)
    ax.scatter(np.array(all_coordinates)[:, 0], np.array(all_coordinates)[:, 1], s=10, c='red')
    ax.scatter(np.array(cluster_coordinates)[:, 0], np.array(cluster_coordinates)[:, 1], s=10, c='blue')
    line_collection = mc.LineCollection(edges, color='orange')
    ax.add_collection(line_collection)

    # TODO adjust axis
    ax.set_ylim([min_ylim - 1, min_ylim + greates_interval + 1])
    ax.set_xlim([min_xlim - 1, min_xlim + greates_interval + 1])
    # ax.set_yticks(np.linspace(min_ylim, max_ylim, num=10))
    # ax.set_xticks(np.linspace(min_xlim, max_xlim, num=10))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # Hide grid lines
    ax.grid(False)

    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])

    ax.invert_xaxis()
    ax.invert_yaxis()

    plt.tight_layout()
    fig.savefig(os.path.join(save_folder, "".join([title, fileformat])))
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
            cyto_color = "#ff7f00"  # orange LICHEN "#ff7f00"
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
        x1, x2 = 0, 1  # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
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
            cyto_color = "#ff7f00"  # orange LICHEN "#ff7f00"
        elif cyto == "IL13":
            cyto_color = "#e41a1c"  # red AE
        else:
            cyto_color = "#377eb8"  # blue PSO
        col_names = "_".join([cyto, 'responder'])

        # Add non-responder and non-cytokine count distribution
        m_cytoresp = (adata.obs["_".join([cyto, 'Responders_counts'])].values > 0) & \
                     (adata.obs["_".join([cyto, 'counts'])].values > 0)
        counts = adata.obs['n_counts'][~m_cytoresp]

        # Apply test to check wether the responders are enriched in the cyto+ spots
        print(corr_stats.apply_wilcoxontest(df_includedcounts[col_names].values, df_excludedcounts[col_names].values))

        # Two-sample Kolmogorov–Smirnov test
        ks_test = scstats.ks_2samp(df_includedcounts[col_names].values, df_excludedcounts[col_names].values)
        print(ks_test)

        included_colname = " ".join([cyto, 'included responder counts'])
        excluded_colname = " ".join([cyto, 'excluded responder counts'])

        df_includedcounts = df_includedcounts.rename(columns={col_names: included_colname})
        df_excludedcounts = df_excludedcounts.rename(columns={col_names: excluded_colname})
        df_otherscounts = pd.DataFrame(counts.values, columns=['Others'])

        df_merged = pd.concat([df_includedcounts[included_colname],
                               df_excludedcounts[excluded_colname], df_otherscounts['Others']], axis=1)
        df_melted = pd.melt(df_merged)

        fig, ax = plt.subplots(figsize=fig_size)
        # Turns off grid on the left Axis.
        ax.grid(False)
        sns.boxplot(x="variable", y="value", data=df_melted, ax=ax, palette=[cyto_color, 'y', 'gray']).set(xlabel='')
        ax.set_ylabel("Counts", fontsize=axis_label_fontsize)

        # statistical annotation
        x1, x2 = 0, 1  # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
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
        plt.savefig(os.path.join(save_folder, "_".join([cyto, 'Boxplots_InExcluded_respondercounts.pdf'])))
        plt.close()


def plot_responder_distributions(adata, df_included_responder, df_excluded_responder, t_cell_cytocines, save_folder):
    """

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


def plot_corrheatmap(df_counts, save_folder):
    """Plot correlation matrix

    Parameters
    ----------
    df_counts : pandas.Dataframe
    save_folder : str

    Returns
    -------

    """
    # remove cluster size, disease and tissue layer
    df_temp = df_counts.drop(['disease', 'Cluster_size', 'tissue_layer'], axis=1, inplace=False)
    # convert all columns of DataFrame
    df_temp = df_temp.apply(pd.to_numeric)

    pcor = pingouin.pairwise_corr(df_temp, covar=df_temp.columns)

    # TODO add permutations to create correlation matrix
    fig, ax = plt.subplots(figsize=fig_size)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    im = ax.matshow(df_temp.corr(), vmin=-1, vmax=1, cmap='coolwarm')
    # im = sns.heatmap(df_temp.corr(), vmin=-1, vmax=1, center=0, ax=ax)

    ax.set_xticks(range(df_temp.shape[1]))
    ax.set_xticklabels(df_temp.columns, rotation=45)
    ax.set_yticks(range(df_temp.shape[1]))
    ax.set_yticklabels(df_temp.columns, rotation=45)

    # Colorbar
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=8)

    plt.tight_layout()
    fig.savefig(os.path.join(save_folder, "_".join(['Correlation_Heatmap', fileformat])))
    plt.close()
