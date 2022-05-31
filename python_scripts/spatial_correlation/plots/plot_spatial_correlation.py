"""Visualise Weighted Spatial Pearson Correlation and Connectivity graphs with and without tissue image
    File name: plot_spatial_correlation.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

# import scripts
from python_scripts.spatial_correlation import corr_statistics as corr_stats
from python_scripts.spatial_correlation import helper_functions, plot_colorbar_legend

# Plotting packages
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns

# System specific
import os

# Calculation packages
import scanpy as sc
import numpy as np
import pandas as pd
import pingouin
from collections import OrderedDict
from scipy import stats


# Figure params
sc.set_figure_params(color_map='viridis')
fig_size, title_fontsize, axis_label_fontsize, legend_fontsize, fileformat, img_key, xy_ticks, text_fontsize = \
    helper_functions.figure_params()


def plot__spatial_correlation(df_counts, cytokine_responders, save_folder, distance):
    """Calculate spatial (Pearson) Correlation between each Cyto+ spot and its nn responder genes spots for each cluster
    -> Plot for Workflow figure (Fig 1)
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
        # Always report the correlation and the p-value:
        # -> Use pearson correlation
        # -> at low statistics the p-value might be infiltrated
        sig_r = pingouin.corr(x=temp_df["_".join([cyto, 'responder'])], y=temp_df[cyto], method='pearson')
        p_vals.append(sig_r['p-val'].values[0])

        #  2. Plot correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        sns.regplot(data=temp_df, x=resp_name, y=cyto, ax=ax, scatter=False, color="black", label=None)
        ax.scatter(data=temp_df, x=resp_name, y=cyto, c='k')

        # Add text: Correlation value and p-value
        ax.text(temp_df.max()[1] / 2 - temp_df.max()[1]/10, temp_df.max()[0],
                'r = {:.2f}; p = {:.2e}'.format(sig_r['r'].values[0], sig_r['p-val'].values[0]),
                fontstyle='italic', fontsize=text_fontsize)

        # Axis params
        ax.set_xlabel(" ".join(["Responder Counts"]), fontsize=axis_label_fontsize)
        if cyto == 'IFNG':
            ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel(" ".join([cyto, 'Counts']), fontsize=axis_label_fontsize)
        ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 50])
        ax.set_ylim([-0.5, temp_df.max()[0] + temp_df.max()[0] / 50])
        ax.tick_params(labelsize=xy_ticks)

        plt.tight_layout()
        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure
        fig.savefig(os.path.join(save_folder, "_".join(['Fig1', str(distance), cyto, resp_name, fileformat])))
        plt.close()

    return p_vals


def plot__swc_tissuelayers(df_counts, cytokine_responders, save_folder, distance):
    """Calculate Correlation between each Cytokine positive spot and its nn responder genes spots for each sample
    -> Plot for Figure 4C

    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    distance : int

    Returns
    -------

    """
    # Get unique combinations of the tissue layer to assign to each a color before the plotting
    # 1. map each string to an integer value
    color_seq = pd.factorize(df_counts['tissue_layer'])[0]
    num_unique_colors = len(np.unique(color_seq))
    # 2. assign a color to each combination
    dict_tissue_color = OrderedDict(zip(df_counts['tissue_layer'], color_seq))
    tissuecomb_colors = plt.cm.get_cmap("tab20", num_unique_colors)

    # 3. Plot Colorbar with tissue layer combinations
    plot_colorbar_legend.plot_standalone_colorbar(
        tissuecomb_colors=tissuecomb_colors, labels=dict_tissue_color.keys(), save_folder=save_folder)

    # 4. get unique cluster sizes
    cluster_sizes = np.unique(df_counts['Cluster_size'])
    plot_colorbar_legend.plot_standalone_legend(cluster_sizes, save_folder=save_folder)

    p_vals = []
    for cyto in cytokine_responders:
        resp_name = "_".join([cyto, 'responder'])
        temp_df = df_counts[[cyto, resp_name, "Cluster_size", "disease", 'tissue_layer']].copy()
        weighted_cytoname = "_".join(['weighted', cyto])

        temp_df = temp_df[~np.isnan(temp_df[cyto].values.astype(np.float64))]
        temp_df[cyto] = temp_df[cyto].values.astype(np.float)
        temp_df[resp_name] = temp_df[resp_name].values.astype(np.float)

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
            x=temp_df[resp_name], y=temp_df[cyto], w=temp_df['Cluster_size'])
        # 1.3 Calculate p-value
        p_w = corr_stats.apply_statstest(df_highcounts=temp_df[resp_name],
                                         df_lowcounts=weighted_cytocounts, correlation_value=correlation)
        p_vals.append(p_w)

        get_sizes = np.unique(temp_df['Cluster_size'])

        # Increase point size
        temp_df["Cluster size"] = (temp_df["Cluster_size"].values**2 * 36).astype(int)

        #  2. Plot correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        ax_sns = sns.regplot(data=temp_df, x=resp_name, y=cyto, ax=ax, scatter=False, color="black", label=None,
                             robust=True)
        plt.setp(ax_sns.lines, zorder=1)
        ax.scatter(data=temp_df, x=resp_name, y=cyto, s="Cluster size", c=cyto_colorseq,
                   cmap=ListedColormap(tissuecomb_colors.colors[np.unique(cyto_colorseq)]), zorder=2)

        # Axis params
        ax.set_xlabel(" ".join(["Responder Counts"]), fontsize=axis_label_fontsize)
        if cyto == 'IFNG':
            ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel(" ".join([cyto, 'Counts']), fontsize=axis_label_fontsize)
        if any(get_sizes > 6):
            ax.set_yticks(np.arange(0, temp_df.max()[0] + 3, 1))
            if temp_df.max()[0] < 10:
                ax.set_yticks(np.arange(0, temp_df.max()[0] + 1, 1))
            else:
                ax.set_yticks(np.arange(0, temp_df.max()[0] + 3, 2))
            ax.set_ylim([-0.5, temp_df.max()[0] + 2])
            ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 13])
            # Add text: Correlation value and p-value
            ax.text(temp_df.max()[1] / 2 - temp_df.max()[1] / 10, temp_df.max()[0] + 1,
                    'r = {:.2f}; p = {:.2e}'.format(correlation, p_w), fontstyle='italic', fontsize=text_fontsize)

        else:
            ax.set_yticks(np.arange(0, temp_df.max()[0] + 1, 1))
            ax.set_ylim([-0.5, temp_df.max()[0] + temp_df.max()[0] / 20])
            ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 20])
            # Add text: Correlation value and p-value
            ax.text(temp_df.max()[1] / 2 - temp_df.max()[1] / 10, temp_df.max()[0],
                    'r = {:.2f}; p = {:.2e}'.format(correlation, p_w), fontstyle='italic', fontsize=text_fontsize)
        ax.tick_params(labelsize=xy_ticks)

        plt.tight_layout()

        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure -> creates figure for Workflow figure 1
        fig.savefig(os.path.join(
            save_folder, "{}_tissuelayers_r{}_{}_{}.{}".format('FigG-I', str(distance), cyto, resp_name, fileformat)),)
        plt.close()

    return p_vals


def plot__swc_disease(df_counts, cytokine_responders, save_folder, distance):
    """Calculate Correlation between each Cytokine positive spot and its nn responder genes spots for each sample
    -> Plot for Figure 4C

    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    distance : int

    Returns
    -------

    """
    # Get unique combinations of the disease to assign to each a color before the plotting
    # 1.1 map each string to an integer value
    color_seq = pd.factorize(df_counts['disease'])[0]
    num_unique_colors = len(np.unique(color_seq))
    # 2.2 assign a color to each combination
    dict_disease_color = OrderedDict(zip(df_counts['disease'], color_seq))
    diseasecomb_colors = plt.cm.get_cmap("tab20", num_unique_colors)

    # 3. Plot Colorbar for diseases
    plot_colorbar_legend.plot_standalone_colorbar_disease(
        diseasecomb_colors=diseasecomb_colors, labels=dict_disease_color.keys(), save_folder=save_folder)

    # 4. get unique cluster sizes
    cluster_sizes = np.unique(df_counts['Cluster_size'])
    plot_colorbar_legend.plot_standalone_legend(cluster_sizes, save_folder=save_folder)

    pvals_corr = {key: [] for key in cytokine_responders.keys()}
    for cyto in cytokine_responders:
        resp_name = "_".join([cyto, 'responder'])
        temp_df = df_counts[[cyto, resp_name, "Cluster_size", "disease", 'tissue_layer']].copy()
        weighted_cytoname = "_".join(['weighted', cyto])

        temp_df = temp_df[~np.isnan(temp_df[cyto].values.astype(np.float64))]
        temp_df[cyto] = temp_df[cyto].values.astype(np.float)
        temp_df[resp_name] = temp_df[resp_name].values.astype(np.float)

        # # map each string to an integer value
        cyto_colorseq = []
        for ind_dict, value in enumerate(temp_df['disease']):
            cyto_colorseq.append(dict_disease_color[value])

        # 1. Preparation phase for weighted correlation
        # Apply weights to cytokine counts and add it to df: Multiply with cluster size
        weighted_cytocounts = temp_df[cyto] * temp_df['Cluster_size']
        temp_df[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)

        #  2. Plot correlation
        get_sizes = []
        pvals_corr[cyto] = {key: [] for key in temp_df['disease'].unique()}

        text_pos = {'Pso': temp_df.max()[1] - temp_df.max()[1] / 8,
                    'AD': temp_df.max()[1] / 2 - temp_df.max()[1] / 14,
                    'LP': temp_df.min()[1]}
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        for ind, diag in enumerate(temp_df['disease'].unique()):
            temp_df_diag = temp_df[temp_df['disease'].str.contains(diag)].copy()
            # Subset color for specific diagnosis
            cyto_colorseq_diag = np.array(cyto_colorseq)[temp_df['disease'].str.contains(diag).values].copy()

            # 2.1 Calculate weighted Pearson Correlation per disease
            if len(temp_df_diag[resp_name]) >= 3:
                correlation = corr_stats.weighted_corr(
                    x=temp_df_diag[resp_name], y=temp_df_diag[cyto], w=temp_df_diag['Cluster_size'])
                # 2.2 Calculate p-value
                p_w = corr_stats.apply_statstest(df_highcounts=temp_df_diag["_".join([cyto, 'responder'])],
                                                 df_lowcounts=weighted_cytocounts, correlation_value=correlation)
            else:
                correlation = np.nan
                p_w = np.nan
            pvals_corr[cyto][diag].append((correlation, p_w))

            get_sizes.extend(np.unique(temp_df_diag['Cluster_size']))

            # Increase point size
            temp_df_diag["Cluster size"] = (temp_df_diag["Cluster_size"].values**2 * 36).astype(int)

            if (len(temp_df_diag[resp_name]) >= 3) & (pvals_corr[cyto][diag][0][1] <= 0.001):
                ax_sns = sns.regplot(
                    data=temp_df_diag, x=resp_name, y=cyto, ax=ax, scatter=False,
                    color=np.unique(ListedColormap(diseasecomb_colors.colors[np.unique(cyto_colorseq_diag)]).colors),
                    label=None, robust=True)
                plt.setp(ax_sns.lines, zorder=1)
                # Add text: Correlation value and p-value per diagnosis
                # x: temp_df.max()[1] / 2 - temp_df.max()[1] / 10, y : temp_df.max()[0] + 1
                ax.text(text_pos[diag], temp_df.max()[0] + 1,
                        r'{}: r = {:.2f}; p = {:.2e}'.format(diag, pvals_corr[cyto][diag][0][0],
                                                             pvals_corr[cyto][diag][0][1]),
                        fontstyle='italic', fontsize=10)

            ax.scatter(data=temp_df_diag, x=resp_name, y=cyto, s="Cluster size", c=cyto_colorseq_diag,
                       cmap=ListedColormap(diseasecomb_colors.colors[np.unique(cyto_colorseq_diag)]), zorder=2,
                       alpha=0.7)

        if any(np.array(get_sizes) > 6):
            ax.set_yticks(np.arange(0, temp_df.max()[0] + 3, 1))
            if temp_df.max()[0] < 10:
                ax.set_yticks(np.arange(0, temp_df.max()[0] + 1, 1))
            else:
                ax.set_yticks(np.arange(0, temp_df.max()[0] + 3, 2))
            ax.set_ylim([-0.5, temp_df.max()[0] + 2])
            ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 13])
        else:
            ax.set_yticks(np.arange(0, temp_df.max()[0] + 1, 1))
            ax.set_ylim([-0.5, temp_df.max()[0] + temp_df.max()[0] / 20])
            ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 20])
        ax.tick_params(labelsize=xy_ticks)

        # Axis params
        ax.set_xlabel(" ".join(["Responder Counts"]), fontsize=axis_label_fontsize)
        if cyto == 'IFNG':
            ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel(" ".join([cyto, 'Counts']), fontsize=axis_label_fontsize)

        plt.tight_layout()

        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure -> creates figure for Workflow figure 1
        fig.savefig(os.path.join(
            save_folder, "{}_diagnosis_r{}_{}_{}.{}".format('FigG-I', str(distance), cyto, resp_name, fileformat)))
        plt.close()

    return pvals_corr


