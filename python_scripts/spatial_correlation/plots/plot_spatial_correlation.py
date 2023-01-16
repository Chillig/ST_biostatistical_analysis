"""Visualise Weighted Spatial Pearson Correlation and Connectivity graphs with and without tissue image
    File name: plot_spatial_correlation.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

# import scripts
from python_scripts.spatial_correlation import corr_statistics as corr_stats
from python_scripts.spatial_correlation import helper_functions
from python_scripts.spatial_correlation.plots import plot_colorbar_legend

# Plotting packages
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import seaborn as sns

# System specific
import os

# Calculation packages
import scanpy as sc
import numpy as np
import pandas as pd
import pingouin
from collections import OrderedDict

# Figure params
fig_size, title_fontsize, axis_label_fontsize, legend_fontsize, fileformat, img_key, xy_ticks, text_fontsize = \
    helper_functions.figure_params()
size_multiplier = 36


def get_correlation_stats(df, resp_name, cyto, weight_name, dict_corr):
    df_wo_nan = df.dropna()

    for method in ['pearson', 'spearman']:
        if method == 'pearson':
            # Weighted Pearson Correlation (WPC)
            correlation = corr_stats.weighted_corr(
                x=df_wo_nan[resp_name], y=df_wo_nan[cyto], w=df_wo_nan[weight_name])
        else:
            # Weighted Spearman Correlation (WSC)
            correlation = corr_stats.calculate_weighted_spearman(
                x=df_wo_nan[resp_name], y=df_wo_nan[cyto], w=df_wo_nan[weight_name])
        # 1.3 Calculate p-value
        p_w = corr_stats.apply_statstest(
            df_highcounts=df_wo_nan[resp_name], df_lowcounts=df_wo_nan[cyto], correlation_value=correlation)

        # Save correlation value and p-value in dictionary
        dict_corr[method].append((cyto, correlation, p_w))

    return dict_corr


def get_disease_correlation_stats(df, resp_name, cyto, weight_name, diag, dict_corr):
    if (len(df[resp_name]) > 2) & (len(df[cyto]) > 2):
        for method in ['pearson', 'spearman']:
            if method == 'pearson':
                # Weighted Pearson Correlation (WPC)
                correlation = corr_stats.weighted_corr(
                    x=df[resp_name], y=df[cyto], w=df[weight_name])
            else:
                # Weighted Spearman Correlation (WSC)
                correlation = corr_stats.calculate_weighted_spearman(
                    x=df[resp_name], y=df[cyto], w=df[weight_name])
            # 1.3 Calculate p-value
            p_w = corr_stats.apply_statstest(
                df_highcounts=df[resp_name], df_lowcounts=df[cyto], correlation_value=correlation)

            # Save correlation value and p-value in dictionary
            dict_corr[method].append((cyto, diag, correlation, p_w))
    else:
        for method in ['pearson', 'spearman']:
            if method == 'pearson':
                # Save correlation value and p-value in dictionary
                dict_corr[method].append((cyto, diag, np.nan, np.nan))
            else:
                # Save correlation value and p-value in dictionary
                dict_corr[method].append((cyto, diag, np.nan, np.nan))

    return dict_corr


def get_linefit(df, resp_name, cyto):
    # Save df in a new df and sort it - otherwise fill_between wont work
    xy_vals = df[[resp_name, cyto]].copy()
    xy_vals = xy_vals.sort_values(by=resp_name)

    # Prepare data for linear regression fit
    x = np.asarray([xy_vals[resp_name].values]).T
    y = np.asarray([xy_vals[cyto].values]).T

    # Fit Line going through origin (0, 0), intercept of 0
    est_fitted, ypred, df_results = corr_stats.create_ls_model(
        x=x, y=y, w=None, const=False)

    return est_fitted, ypred, df_results, x


def plot__stwc_tissuelayers(df_counts, cytokine_responders, save_folder, distance, corr_method):
    """Calculate Correlation between each Cytokine positive spot and its nn responder genes spots for each sample
    -> Plot for Figure 4G-I

    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    distance : int
    corr_method : str

    Returns
    -------

    """
    # Get unique combinations of the tissue layer to assign to each a color before the plotting
    # 1. map each string to an integer value
    # color_seq = pd.factorize(df_counts['tissue_layer'])[0]
    # num_unique_colors = len(np.unique(color_seq))
    # 2. assign a color to each combination
    dict_tissue_color = {'upper EPIDERMIS': 0,
                         'basal EPIDERMIS': 1,
                         'middle EPIDERMIS': 2,
                         'middle EPIDERMIS_upper EPIDERMIS': 3,
                         'basal EPIDERMIS_upper EPIDERMIS': 4,
                         'basal EPIDERMIS_middle EPIDERMIS': 5,
                         'basal EPIDERMIS_middle EPIDERMIS_upper EPIDERMIS': 6}
    num_unique_colors = len(dict_tissue_color.values())
    tissuecomb_colors = plt.cm.get_cmap("tab20", num_unique_colors)

    # 3. Plot Colorbar with tissue layer combinations
    plot_colorbar_legend.plot_standalone_colorbar(
        tissuecomb_colors=tissuecomb_colors, labels=dict_tissue_color.keys(), save_folder=save_folder)

    # 4. get unique cluster sizes
    cluster_sizes = np.unique(df_counts['Cluster_size'])
    plot_colorbar_legend.plot_standalone_legend(cluster_sizes, save_folder=save_folder, size_multiplier=size_multiplier)

    dict_corr = dict({'pearson': [], 'spearman': []})
    for cyto in cytokine_responders:
        resp_name = "_".join([cyto, 'responder'])
        temp_df = df_counts[[cyto, resp_name, "Cluster_size", "disease", 'tissue_layer', 'Cluster_num_spots']].copy()
        weighted_cytoname = "_".join(['weighted', cyto])

        temp_df = temp_df[~np.isnan(temp_df[cyto].values.astype(np.float64))]
        temp_df[cyto] = temp_df[cyto].values.astype(np.float)
        temp_df[resp_name] = temp_df[resp_name].values.astype(np.float)

        # # map each string to an integer value
        cyto_colorseq = []
        for ind_dict, value in enumerate(temp_df['tissue_layer']):
            cyto_colorseq.append(dict_tissue_color[value])

        # 1. Calculate correlation
        # 1.1 Apply weights to cytokine counts and add it to df
        weighted_cytocounts = temp_df[cyto]
        temp_df[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)

        dict_corr = get_correlation_stats(
            df=temp_df, resp_name=resp_name, cyto=cyto, weight_name=weighted_cytoname, dict_corr=dict_corr)

        get_sizes = np.unique(temp_df['Cluster_size'])
        if any(get_sizes > 6):
            increased_ylimit = 3
        else:
            increased_ylimit = 1

        # Increase point size
        temp_df["Cluster size"] = (temp_df["Cluster_size"].values**2 * size_multiplier).astype(int)

        # Fit Line going through origin (0, 0), intercept of 0
        est_woweights_fitted, ypred_woweights, df_woweights, x = get_linefit(
            df=temp_df, resp_name=resp_name, cyto=cyto)

        #  Save in .csv
        temp_df['weights'] = temp_df[cyto]

        # 2. Plot Correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        # Unweighted: Plot weighted line fit
        ax.plot(x, ypred_woweights, label="fit", color='black', zorder=1, lw=2)
        # Unweighted: plotting error band
        ax.fill_between(x.transpose()[0], df_woweights['mean_ci_lower'].values, df_woweights['mean_ci_upper'].values,
                        alpha=.1, label='5 - sigma interval', color='black', lw=0.1)
        # Plot unweighted points
        ax.scatter(data=temp_df, x=resp_name, y=cyto, c=cyto_colorseq, s='Cluster size', alpha=0.8,
                   cmap=ListedColormap(tissuecomb_colors.colors[np.unique(cyto_colorseq)]), zorder=2)

        # Axis params
        xlabel_name = resp_name.split('_')
        if len(xlabel_name) > 2:
            if xlabel_name[1] == 'IFNG':
                ax.set_xlabel(r'IFN-$\gamma$ Responder Counts', fontsize=axis_label_fontsize)
            else:
                ax.set_xlabel("{} Responder Counts".format(xlabel_name[1]), fontsize=axis_label_fontsize)
        else:
            ax.set_xlabel("Responder Counts", fontsize=axis_label_fontsize)
        if cyto.split('_')[0] == 'IFNG':
            ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel(" ".join([cyto.split('_')[0], 'Counts']), fontsize=axis_label_fontsize)

        # Add correlation value and p-value
        corr_pval = [a_tuple for a_tuple in dict_corr[corr_method] if a_tuple[0] == cyto]

        if temp_df[cyto].max() < 10:
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 1))
        elif (temp_df[cyto].max() >= 10) & (temp_df[cyto].max() < 50):
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 5))
        elif (temp_df[cyto].max() >= 50) & (temp_df[cyto].max() < 100):
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 10))
        elif (temp_df[cyto].max() >= 100) & (temp_df[cyto].max() < 1000):
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 100))
        else:
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 20))

        if any(get_sizes > 6):
            ax.set_ylim([-0.1, temp_df[cyto].max() + 2])
            ax.set_xlim([0.5, temp_df[resp_name].max() + temp_df[resp_name].max() / 13])
            # Add text: Correlation value and p-value
            ax.text(temp_df[resp_name].max() / 2 - temp_df[resp_name].max() / 10,
                    temp_df[cyto].max() + 1,
                    'r = {:.2f}; p = {:.2e}'.format(corr_pval[0][1], corr_pval[0][2]),
                    fontstyle='italic', fontsize=text_fontsize)

        else:
            ax.set_ylim([-0.1, temp_df[cyto].max() + temp_df[cyto].max() / 20])
            ax.set_xlim([0.5, temp_df[resp_name].max() + temp_df[resp_name].max() / 20])
            # Add text: Correlation value and p-value
            ax.text(temp_df[resp_name].max() / 2 -
                    temp_df[resp_name].max() / 10,
                    temp_df[cyto].max(),
                    'r = {:.2f}; p = {:.2e}'.format(corr_pval[0][1], corr_pval[0][2]),
                    fontstyle='italic', fontsize=text_fontsize)
        ax.tick_params(labelsize=xy_ticks)

        plt.tight_layout()

        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure -> creates figure for Workflow figure 1
        fig.savefig(os.path.join(
            save_folder, "{}_tissuelayers_weighted_transcripts_{}_r{}_{}_{}{}".format(
                'FigF-H', corr_method, str(distance), cyto, resp_name, fileformat)),)
        plt.close()

    return dict_corr


def plot__stwc_patients(df_counts, cytokine_responders, save_folder, distance, corr_method):
    """Calculate Correlation between each Cytokine positive spot and its nn responder genes spots for each sample
    -> Plot for Figure 4G-I

    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    distance : int
    corr_method : str

    Returns
    -------

    """
    # Get unique combinations of the tissue layer to assign to each a color before the plotting
    # 1. assign a color to each combination
    color_seq = pd.factorize(df_counts['Patient'].astype(int))[0]
    num_unique_colors = len(np.unique(color_seq))
    # 1.1 assign a color to each combination
    dict_patient_color = OrderedDict(zip(df_counts['Patient'].astype(int), color_seq))
    patient_colors = mpl.colors.ListedColormap(sc.pl.palettes.default_102[:num_unique_colors])

    # 2. Plot Colorbar with tissue layer combinations
    plot_colorbar_legend.plot_standalone_colorbar_patient(
        patient_colors=patient_colors, labels=dict_patient_color.keys(), save_folder=save_folder)

    # 3. get unique cluster sizes
    cluster_sizes = np.unique(df_counts['Cluster_size'])
    plot_colorbar_legend.plot_standalone_legend(cluster_sizes, save_folder=save_folder, size_multiplier=size_multiplier)

    dict_corr = dict({'pearson': [], 'spearman': []})
    for cyto in cytokine_responders:
        resp_name = "_".join([cyto, 'responder'])
        temp_df = df_counts[[
            cyto, resp_name, "Cluster_size", "disease", 'tissue_layer', 'Patient', 'Cluster_num_spots']].copy()
        weighted_cytoname = "_".join(['weighted', cyto])

        temp_df = temp_df[~np.isnan(temp_df[cyto].values.astype(np.float64))]
        temp_df[cyto] = temp_df[cyto].values.astype(np.float)
        temp_df[resp_name] = temp_df[resp_name].values.astype(np.float)

        # # map each string to an integer value
        cyto_colorseq = []
        for ind_dict, value in enumerate(temp_df['Patient']):
            cyto_colorseq.append(dict_patient_color[value])

        # 1. Calculate correlation
        # 1.1 Apply weights to cytokine counts and add it to df
        weighted_cytocounts = temp_df[cyto]
        temp_df[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)

        dict_corr = get_correlation_stats(
            df=temp_df, resp_name=resp_name, cyto=cyto, weight_name=weighted_cytoname, dict_corr=dict_corr)

        get_sizes = np.unique(temp_df['Cluster_size'])
        if any(get_sizes > 6):
            increased_ylimit = 3
        else:
            increased_ylimit = 1

        # Increase point size
        temp_df["Cluster size"] = (temp_df["Cluster_size"].values**2 * size_multiplier).astype(int)

        # Fit Line going through origin (0, 0), intercept of 0
        est_woweights_fitted, ypred_woweights, df_woweights, x = get_linefit(
            df=temp_df, resp_name=resp_name, cyto=cyto)

        #  Save in .csv
        temp_df['weights'] = temp_df[cyto]

        # 2. Plot Correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        # Unweighted: Plot weighted line fit
        ax.plot(x, ypred_woweights, label="fit", color='black', zorder=1, lw=2)
        # Unweighted: plotting error band
        ax.fill_between(x.transpose()[0], df_woweights['mean_ci_lower'].values, df_woweights['mean_ci_upper'].values,
                        alpha=.1, label='5 - sigma interval', color='black', lw=0.1)

        # Plot unweighted points
        ax.scatter(data=temp_df, x=resp_name, y=cyto, c=cyto_colorseq, s='Cluster size', alpha=0.8,
                   cmap=ListedColormap(np.array(patient_colors.colors)[np.unique(cyto_colorseq)]), zorder=2)

        # Axis params
        xlabel_name = resp_name.split('_')
        if len(xlabel_name) > 2:
            if xlabel_name[1] == 'IFNG':
                ax.set_xlabel(r'IFN-$\gamma$ Responder Counts', fontsize=axis_label_fontsize)
            else:
                ax.set_xlabel("{} Responder Counts".format(xlabel_name[1]), fontsize=axis_label_fontsize)
        else:
            ax.set_xlabel("Responder Counts", fontsize=axis_label_fontsize)
        if cyto.split('_')[0] == 'IFNG':
            ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel(" ".join([cyto.split('_')[0], 'Counts']), fontsize=axis_label_fontsize)

        # Add correlation value and p-value
        corr_pval = [a_tuple for a_tuple in dict_corr[corr_method] if a_tuple[0] == cyto]

        if temp_df[cyto].max() < 10:
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 1))
        elif (temp_df[cyto].max() >= 10) & (temp_df[cyto].max() < 50):
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 5))
        elif (temp_df[cyto].max() >= 50) & (temp_df[cyto].max() < 100):
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 10))
        elif (temp_df[cyto].max() >= 100) & (temp_df[cyto].max() < 1000):
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 100))
        else:
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 20))

        if any(get_sizes > 6):
            ax.set_ylim([-0.1, temp_df[cyto].max() + 2])
            ax.set_xlim([0.5, temp_df[resp_name].max() + temp_df[resp_name].max() / 13])
            # Add text: Correlation value and p-value
            ax.text(temp_df[resp_name].max() / 2 - temp_df[resp_name].max() / 10,
                    temp_df[cyto].max() + 1,
                    'r = {:.2f}; p = {:.2e}'.format(corr_pval[0][1], corr_pval[0][2]),
                    fontstyle='italic', fontsize=text_fontsize)

        else:
            ax.set_ylim([-0.1, temp_df[cyto].max() + temp_df[cyto].max() / 20])
            ax.set_xlim([0.5, temp_df[resp_name].max() + temp_df[resp_name].max() / 20])
            # Add text: Correlation value and p-value
            ax.text(temp_df[resp_name].max() / 2 -
                    temp_df[resp_name].max() / 10,
                    temp_df[cyto].max(),
                    'r = {:.2f}; p = {:.2e}'.format(corr_pval[0][1], corr_pval[0][2]),
                    fontstyle='italic', fontsize=text_fontsize)
        ax.tick_params(labelsize=xy_ticks)

        plt.tight_layout()

        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure -> creates figure for Workflow figure 1
        fig.savefig(os.path.join(
            save_folder, "{}_patients_weighted_transcripts_{}_r{}_{}_{}{}".format(
                'FigS10', corr_method, str(distance), cyto, resp_name, fileformat)),)
        plt.close()

    return dict_corr


def plot__stwc_tissuelayers_figure1g(df_counts, cytokine_responders, save_folder, distance, corr_method):
    """Calculate Correlation between each Cytokine positive spot and its nn responder genes spots for each sample
    -> Plot for Figure 4G-I

    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders :
    save_folder : str
    distance : int
    corr_method : str

    Returns
    -------

    """
    # Get unique combinations of the tissue layer to assign to each a color before the plotting
    # 1. assign a color to each combination
    dict_tissue_color = {'upper EPIDERMIS': 0,
                         'basal EPIDERMIS': 1,
                         'middle EPIDERMIS': 2,
                         'middle EPIDERMIS_upper EPIDERMIS': 3,
                         'basal EPIDERMIS_upper EPIDERMIS': 4,
                         'basal EPIDERMIS_middle EPIDERMIS': 5,
                         'basal EPIDERMIS_middle EPIDERMIS_upper EPIDERMIS': 6}
    num_unique_colors = len(dict_tissue_color.values())
    tissuecomb_colors = plt.cm.get_cmap("tab20", num_unique_colors)

    # 2. Plot Colorbar with tissue layer combinations
    plot_colorbar_legend.plot_standalone_colorbar(
        tissuecomb_colors=tissuecomb_colors, labels=dict_tissue_color.keys(), save_folder=save_folder)

    dict_corr = dict({'pearson': [], 'spearman': []})
    for cyto in cytokine_responders:
        resp_name = "_".join([cyto, 'responder'])
        temp_df = df_counts[[cyto, resp_name, "Cluster_size", "disease", 'tissue_layer', 'Cluster_num_spots']].copy()
        weighted_cytoname = "_".join(['weighted', cyto])

        temp_df = temp_df[~np.isnan(temp_df[cyto].values.astype(np.float64))]
        temp_df[cyto] = temp_df[cyto].values.astype(np.float)
        temp_df[resp_name] = temp_df[resp_name].values.astype(np.float)

        # # map each string to an integer value
        cyto_colorseq = []
        for ind_dict, value in enumerate(temp_df['tissue_layer']):
            cyto_colorseq.append(dict_tissue_color[value])

        # 1. Calculate correlation
        # 1.1 Apply weights to cytokine counts and add it to df
        weighted_cytocounts = temp_df[cyto]
        temp_df[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)

        get_sizes = np.unique(temp_df['Cluster_size'])
        if any(get_sizes > 6):
            increased_ylimit = 3
        else:
            increased_ylimit = 1

        # Fit Line going through origin (0, 0), intercept of 0
        est_woweights_fitted, ypred_woweights, df_woweights, x = get_linefit(
            df=temp_df, resp_name=resp_name, cyto=cyto)

        # 1. Calculate correlation
        # 1.1 Apply weights to cytokine counts and add it to df
        temp_df['weights'] = temp_df[cyto]
        temp_df[weighted_cytoname] = temp_df['weights'].values.astype(np.float64)

        dict_corr = get_correlation_stats(
            df=temp_df, resp_name=resp_name, cyto=cyto, weight_name=weighted_cytoname, dict_corr=dict_corr)

        # Add correlation value and p-value
        corr_pval = [a_tuple for a_tuple in dict_corr[corr_method] if a_tuple[0] == cyto]

        # 2. Plot Correlation
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.grid(False)
        # Unweighted: Plot weighted line fit
        ax.plot(x, ypred_woweights, label="fit", color='black', zorder=1, lw=2)
        # Unweighted: plotting error band
        ax.fill_between(x.transpose()[0], df_woweights['mean_ci_lower'].values, df_woweights['mean_ci_upper'].values,
                        alpha=.1, label='5 - sigma interval', color='black', lw=0.1)
        # Plot unweighted points
        ax.scatter(data=temp_df, x=resp_name, y=cyto, c='black', alpha=1, s=100,
                   cmap=ListedColormap(tissuecomb_colors.colors[np.unique(cyto_colorseq)]),
                   zorder=2, edgecolors='black')

        # Axis params
        xlabel_name = resp_name.split('_')
        if len(xlabel_name) > 2:
            ax.set_xlabel("{} Responder Counts".format(xlabel_name[1]), fontsize=18)
        else:
            ax.set_xlabel("Responder Counts", fontsize=18)
        ax.set_ylabel(" ".join([cyto.split('_')[0], 'Counts']), fontsize=18)

        if temp_df[cyto].max() < 10:
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 1))
        elif (temp_df[cyto].max() >= 10) & (temp_df[cyto].max() < 20):
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 2))
        elif (temp_df[cyto].max() >= 20) & (temp_df[cyto].max() < 30):
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 3))
        elif (temp_df[cyto].max() >= 30) & (temp_df[cyto].max() < 50):
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 5))
        elif (temp_df[cyto].max() >= 50) & (temp_df[cyto].max() < 100):
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 10))
        elif (temp_df[cyto].max() >= 100) & (temp_df[cyto].max() < 1000):
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 100))
        else:
            ax.set_yticks(np.arange(0, temp_df[cyto].max() + increased_ylimit, 20))

        # Add text: Correlation value and p-value
        ax.text(temp_df[resp_name].max() / 2 - temp_df[resp_name].max() / 5,
                temp_df[cyto].max() + 1,
                'r = {:.2f}; p = {:.2e}'.format(corr_pval[0][1], corr_pval[0][2]),
                fontstyle='italic', fontsize=text_fontsize)

        ax.tick_params(labelsize=16)

        plt.tight_layout()

        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure -> creates figure for Workflow figure 1
        fig.savefig(os.path.join(
            save_folder, "{}_tissuelayers_weighted_transcripts_{}_r{}_{}_{}{}".format(
                'Fig1G', corr_method, str(distance), cyto, resp_name, fileformat)),)
        plt.close()
