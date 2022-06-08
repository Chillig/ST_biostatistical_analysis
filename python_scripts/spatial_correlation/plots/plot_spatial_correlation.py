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


def plot__spatial_correlation(df_counts, cytokine_responders, save_folder, distance, corr_method: str):
    """Calculate spatial (Pearson) Correlation between each Cyto+ spot and its nn responder genes spots for each cluster
    -> Plot for Workflow figure (Fig 1)
    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    distance : int
    corr_method : str

    Returns
    -------
    list of p-values

    """

    dict_corr = dict({'pearson': [], 'spearman': []})
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
        for method in ['pearson', 'spearman']:
            sig_r = pingouin.corr(x=temp_df[resp_name], y=temp_df[cyto], method=method)
            dict_corr[method].append((cyto, sig_r['r'].values[0], sig_r['p-val'].values[0]))

        #  2. Plot correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        sns.regplot(data=temp_df, x=resp_name, y=cyto, ax=ax, scatter=False, color="black", label=None)
        ax.scatter(data=temp_df, x=resp_name, y=cyto, c='k')

        # Add text: Correlation value and p-value
        corr_pval = [a_tuple for a_tuple in dict_corr[corr_method] if a_tuple[0] == cyto]

        ax.text(temp_df.max()[1] / 2 - temp_df.max()[1]/10, temp_df.max()[0],
                'r = {:.2f}; p = {:.2e}'.format(corr_pval[0][1], corr_pval[0][2]),
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

    return dict_corr


def plot__spatial_correlation__wdisease(df_counts, cytokine_responders, save_folder, distance, corr_method: str):
    """Calculate spatial (Pearson) Correlation between each Cyto+ spot and its nn responder genes spots for each cluster
    -> Plot for Workflow figure (Fig 1)
    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    distance : int
    corr_method : str

    Returns
    -------
    list of p-values

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

    dict_corr = dict({'pearson': [], 'spearman': []})
    for cyto in cytokine_responders:
        resp_name = "_".join([cyto, 'responder'])
        weighted_cytoname = "_".join(['weighted', cyto])
        temp_df = df_counts[[cyto, resp_name, "Cluster_size", "disease", 'tissue_layer']].copy()

        temp_df = temp_df[~np.isnan(temp_df[cyto].values.astype(np.float64))]
        temp_df[cyto] = temp_df[cyto].values.astype(np.float)
        temp_df[resp_name] = temp_df[resp_name].values.astype(np.float)

        # # map each string to an integer value
        cyto_colorseq = []
        for ind_dict, value in enumerate(temp_df['disease']):
            cyto_colorseq.append(dict_disease_color[value])

        # 1. Plot correlation
        text_pos = {'Pso': temp_df.max()[1] - temp_df.max()[1] / 8,
                    'AD': temp_df.max()[1] / 2 - temp_df.max()[1] / 14,
                    'LP': temp_df.min()[1]}

        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        for ind, diag in enumerate(temp_df['disease'].unique()):
            temp_df_diag = temp_df[temp_df['disease'].str.contains(diag)].copy()
            # Subset color for specific diagnosis
            cyto_colorseq_diag = np.array(cyto_colorseq)[temp_df['disease'].str.contains(diag).values].copy()

            # 1. Apply weights to cytokine counts and add it to df
            weighted_cytocounts = temp_df_diag[cyto]
            temp_df_diag[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)

            # 2. Calculate Correlations
            dict_corr = get_disease_correlation_stats(
                df=temp_df_diag, resp_name=resp_name, cyto=cyto, weight_name=weighted_cytoname, dict_corr=dict_corr,
                diag=diag)

            # Read out correlation value and p-value for each disease
            corr_pval = [a_tuple for a_tuple in dict_corr[corr_method]
                         if (a_tuple[0] == cyto) and (a_tuple[1] == diag)]

            # 3. Fit Line going through origin (0, 0), intercept of 0
            est_woweights_fitted, ypred_woweights, df_woweights, x = get_linefit(
                df=temp_df_diag, resp_name=resp_name, cyto=cyto)

            # Add Correlation value
            if (len(temp_df_diag[resp_name]) >= 3) & (corr_pval[0][3] <= 0.01):
                # Unweighted: Plot weighted line fit
                ax.plot(
                    x, ypred_woweights, label="fit",
                    color=np.unique(ListedColormap(diseasecomb_colors.colors[np.unique(cyto_colorseq_diag)]).colors),
                    zorder=1, lw=2)
                # Unweighted: plotting error band
                ax.fill_between(
                    x.transpose()[0], df_woweights['mean_ci_lower'].values, df_woweights['mean_ci_upper'].values,
                    alpha=.1, label='5 - sigma interval', lw=0.1,
                    color=np.unique(ListedColormap(diseasecomb_colors.colors[np.unique(cyto_colorseq_diag)]).colors))
                # Add text: Correlation value and p-value per diagnosis
                # x: temp_df.max()[1] / 2 - temp_df.max()[1] / 10, y : temp_df.max()[0] + 1
                ax.text(text_pos[diag], temp_df.max()[0] + 1,
                        r'{}: r = {:.2f}; p = {:.2e}'.format(diag, corr_pval[0][2], corr_pval[0][3]),
                        fontstyle='italic', fontsize=10)

            ax.scatter(data=temp_df_diag, x=resp_name, y=cyto, c=cyto_colorseq_diag,
                       cmap=ListedColormap(diseasecomb_colors.colors[np.unique(cyto_colorseq_diag)]), zorder=2,
                       alpha=0.7)

        # Axis params
        ax.set_xlabel(" ".join(["Responder Counts"]), fontsize=axis_label_fontsize)
        if cyto == 'IFNG':
            ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel(" ".join([cyto, 'Counts']), fontsize=axis_label_fontsize)

        plt.tight_layout()

        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        ax.set_ylim([-0.1, temp_df[cyto].max() + 2])
        ax.set_xlim([0.5, temp_df[resp_name].max() + temp_df[resp_name].max() / 13])
        ax.tick_params(labelsize=xy_ticks)

        # 3. Save figure -> creates figure for Workflow figure 1
        fig.savefig(os.path.join(
            save_folder, "{}_diagnosis_{}_r{}_{}_{}.{}".format(
                'Fig1', corr_method, str(distance), cyto, resp_name, fileformat)))
        plt.close()

    return dict_corr


def plot__swc_tissuelayers(df_counts, cytokine_responders, save_folder, distance, corr_method):
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
    plot_colorbar_legend.plot_standalone_legend(cluster_sizes, save_folder=save_folder, size_multiplier=size_multiplier)

    dict_corr = dict({'pearson': [], 'spearman': []})
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
        # 1.2 Calculate Correlations
        for method in ['pearson', 'spearman']:
            if method == 'pearson':
                # Weighted Pearson Correlation (WPC)
                correlation = corr_stats.weighted_corr(
                    x=temp_df[resp_name], y=temp_df[cyto], w=temp_df['Cluster_size'])
            else:
                # Weighted Spearman Correlation (WSC)
                correlation = corr_stats.calculate_weighted_spearman(
                    x=temp_df[resp_name], y=temp_df[cyto], w=temp_df['Cluster_size'])
            # 1.3 Calculate p-value
            p_w = corr_stats.apply_statstest(df_highcounts=temp_df[resp_name],
                                             df_lowcounts=weighted_cytocounts, correlation_value=correlation)
            dict_corr[method].append((cyto, correlation, p_w))

        get_sizes = np.unique(temp_df['Cluster_size'])

        # Increase point size
        temp_df["Cluster size"] = (temp_df["Cluster_size"].values**2 * size_multiplier).astype(int)
        temp_df[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)

        #  2. Plot correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        # TODO adjust line fitting
        ax_sns = sns.regplot(data=temp_df, x=resp_name, y=cyto, ax=ax, scatter=False, color="black", label=None,
                             robust=True)
        plt.setp(ax_sns.lines, zorder=1)
        # TODO plot weighted points
        ax.scatter(data=temp_df, x=resp_name, y=cyto, s="Cluster size", c=cyto_colorseq,
                   cmap=ListedColormap(tissuecomb_colors.colors[np.unique(cyto_colorseq)]), zorder=2)

        # Axis params
        ax.set_xlabel(" ".join(["Responder Counts"]), fontsize=axis_label_fontsize)
        if cyto == 'IFNG':
            ax.set_ylabel(r'IFN-$\gamma$ Counts', fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel(" ".join([cyto, 'Counts']), fontsize=axis_label_fontsize)

        # Add correlation value and p-value
        corr_pval = [a_tuple for a_tuple in dict_corr[corr_method] if a_tuple[0] == cyto]

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
                    'r = {:.2f}; p = {:.2e}'.format(corr_pval[0][1], corr_pval[0][2]),
                    fontstyle='italic', fontsize=text_fontsize)

        else:
            ax.set_yticks(np.arange(0, temp_df.max()[0] + 1, 1))
            ax.set_ylim([-0.5, temp_df.max()[0] + temp_df.max()[0] / 20])
            ax.set_xlim([-0.5, temp_df.max()[1] + temp_df.max()[1] / 20])
            # Add text: Correlation value and p-value
            ax.text(temp_df.max()[1] / 2 - temp_df.max()[1] / 10, temp_df.max()[0],
                    'r = {:.2f}; p = {:.2e}'.format(corr_pval[0][1], corr_pval[0][2]),
                    fontstyle='italic', fontsize=text_fontsize)
        ax.tick_params(labelsize=xy_ticks)

        plt.tight_layout()

        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure -> creates figure for Workflow figure 1
        fig.savefig(os.path.join(
            save_folder, "{}_tissuelayers_{}_r{}_{}_{}.{}".format(
                'FigG-I', corr_method, str(distance), cyto, resp_name, fileformat)),)
        plt.close()

    return dict_corr


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
        # temp_df[[cyto, resp_name, 'weights', "Cluster_size", "disease", 'tissue_layer']].to_csv(
        #     os.path.join(save_folder, 'radius{}_cytokine{}_vs_{}.csv'.format(
        #         str(distance), cyto, resp_name)))

        # 2. Plot Correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        # Unweighted: Plot weighted line fit
        ax.plot(x, ypred_woweights, label="fit", color='black', zorder=1, lw=2)
        # Unweighted: plotting error band
        ax.fill_between(x.transpose()[0], df_woweights['mean_ci_lower'].values, df_woweights['mean_ci_upper'].values,
                        alpha=.1, label='5 - sigma interval', color='black', lw=0.1)
        # TODO add dashed fitted line of bulk RNA-seq data
        # # Weighted: Plot weighted line fit
        # ax.plot(x, ypred_wweights, label="fit", color='black', zorder=1, lw=2)
        # # weighted: plotting error band
        # ax.fill_between(x.transpose()[0], df_wweights['mean_ci_lower'].values, df_wweights['mean_ci_upper'].values,
        #                 alpha=.1, label='5 - sigma interval', color='black', lw=0.1)
        # weighted: plotting the confidence intervals of 1 sigma
        # ax.fill_between(xy_vals[resp_name], bound_lower, bound_upper, color='orange', alpha=0.1)
        # Plot unweighted points
        ax.scatter(data=temp_df, x=resp_name, y=cyto, c=cyto_colorseq, s='Cluster size',
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

        # # Add text: slope and intercept term
        # ax.text(temp_df[resp_name].max() / 2 + temp_df[resp_name].max() / 10, 0,
        #         'm = {:.2e}; b = {:.2f}'.format(est_woweights_fitted.params[0], est_woweights_fitted.params[1]),
        #         fontstyle='italic', fontsize=text_fontsize)

        plt.tight_layout()

        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure -> creates figure for Workflow figure 1
        fig.savefig(os.path.join(
            save_folder, "{}_tissuelayers_weighted_transcripts_{}_r{}_{}_{}{}".format(
                'FigF-H', corr_method, str(distance), cyto, resp_name, fileformat)),)
        plt.close()

    return dict_corr


def plot__swc_disease(df_counts, cytokine_responders, save_folder, corr_method, distance):
    """Calculate Correlation between each Cytokine positive spot and its nn responder genes spots for each sample
    -> Plot for Figure 4G-I

    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    corr_method: str
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
    plot_colorbar_legend.plot_standalone_legend(cluster_sizes, save_folder=save_folder, size_multiplier=size_multiplier)

    dict_corr = dict({'pearson': [], 'spearman': []})
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
            # 1.2 Calculate Correlations
            for method in ['pearson', 'spearman']:
                if len(temp_df_diag[resp_name]) >= 3:
                    if method == 'pearson':
                        # Weighted Pearson Correlation (WPC)
                        correlation = corr_stats.weighted_corr(
                            x=temp_df_diag[resp_name], y=temp_df_diag[cyto], w=temp_df_diag['Cluster_size'])
                    else:
                        # Weighted Spearman Correlation (WSC)
                        correlation = corr_stats.calculate_weighted_spearman(
                            x=temp_df_diag[resp_name], y=temp_df_diag[cyto], w=temp_df_diag['Cluster_size'])
                    # 1.3 Calculate p-value
                    p_w = corr_stats.apply_statstest(df_highcounts=temp_df_diag[resp_name],
                                                     df_lowcounts=weighted_cytocounts, correlation_value=correlation)
                else:
                    correlation = np.nan
                    p_w = np.nan
                dict_corr[method].append((cyto, diag, correlation, p_w))

            get_sizes.extend(np.unique(temp_df_diag['Cluster_size']))

            # Increase point size
            temp_df_diag["Cluster size"] = (temp_df_diag["Cluster_size"].values**2 * size_multiplier).astype(int)

            # Read out correlation value and p-value for each disease
            corr_pval = [a_tuple for a_tuple in dict_corr[corr_method]
                         if (a_tuple[0] == cyto) and (a_tuple[1] == diag)]

            # Add Correlation value
            if (len(temp_df_diag[resp_name]) >= 3) & (corr_pval[0][3] <= 0.01):
                ax_sns = sns.regplot(
                    data=temp_df_diag, x=resp_name, y=cyto, ax=ax, scatter=False,
                    color=np.unique(ListedColormap(diseasecomb_colors.colors[np.unique(cyto_colorseq_diag)]).colors),
                    label=None, robust=True)
                plt.setp(ax_sns.lines, zorder=1)
                # Add text: Correlation value and p-value per diagnosis
                # x: temp_df.max()[1] / 2 - temp_df.max()[1] / 10, y : temp_df.max()[0] + 1
                ax.text(text_pos[diag], temp_df.max()[0] + 1,
                        r'{}: r = {:.2f}; p = {:.2e}'.format(diag, corr_pval[0][2], corr_pval[0][3]),
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
            save_folder, "{}_diagnosis_{}_r{}_{}_{}.{}".format(
                'FigG-I', corr_method, str(distance), cyto, resp_name, fileformat)))
        plt.close()

    return dict_corr


def plot__swtc_disease(df_counts, cytokine_responders, save_folder, corr_method, distance):
    """Calculate weighted by transcripts Correlation between each Cytokine positive spot and its
    nn responder genes spots for each sample
    -> Plot for Figure 4G-I

    Parameters
    ----------
    df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    corr_method: str
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
    plot_colorbar_legend.plot_standalone_legend(cluster_sizes, save_folder=save_folder, size_multiplier=size_multiplier)

    dict_corr = dict({'pearson': [], 'spearman': []})
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

        #  1. Plot correlation
        get_sizes = []
        text_pos = {'Pso': temp_df.max()[1] - temp_df.max()[1] / 8,
                    'AD': temp_df.max()[1] / 2 - temp_df.max()[1] / 14,
                    'LP': temp_df.min()[1]}

        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        for ind, diag in enumerate(temp_df['disease'].unique()):
            temp_df_diag = temp_df[temp_df['disease'].str.contains(diag)].copy()
            # Subset color for specific diagnosis
            cyto_colorseq_diag = np.array(cyto_colorseq)[temp_df['disease'].str.contains(diag).values].copy()

            # 2.1 Preparation phase for weighted correlation
            # Apply weights to cytokine counts and add it to df: Multiply with cluster size
            weighted_cytocounts = temp_df_diag[cyto] * temp_df_diag[cyto]
            temp_df_diag[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)

            # 2.1 Calculate weighted Pearson Correlation per disease
            # Calculate Correlations
            for method in ['pearson', 'spearman']:
                if len(temp_df_diag[resp_name]) >= 3:
                    if method == 'pearson':
                        # Weighted Pearson Correlation (WPC)
                        correlation = corr_stats.weighted_corr(
                            x=temp_df_diag[resp_name], y=temp_df_diag[cyto], w=temp_df_diag[cyto])
                    else:
                        # Weighted Spearman Correlation (WSC)
                        correlation = corr_stats.calculate_weighted_spearman(
                            x=temp_df_diag[resp_name], y=temp_df_diag[cyto], w=temp_df_diag[cyto])
                    # 1.3 Calculate p-value
                    p_w = corr_stats.apply_statstest(df_highcounts=temp_df_diag[resp_name],
                                                     df_lowcounts=weighted_cytocounts, correlation_value=correlation)
                else:
                    correlation = np.nan
                    p_w = np.nan
                dict_corr[method].append((cyto, diag, correlation, p_w))

            get_sizes.extend(np.unique(temp_df_diag['Cluster_size']))

            # Increase point size
            temp_df_diag["Cluster size"] = (temp_df_diag["Cluster_size"].values**2 * size_multiplier).astype(int)

            # Read out correlation value and p-value for each disease
            corr_pval = [a_tuple for a_tuple in dict_corr[corr_method]
                         if (a_tuple[0] == cyto) and (a_tuple[1] == diag)]

            # Add Correlation value
            if (len(temp_df_diag[resp_name]) >= 3) & (corr_pval[0][3] <= 0.01):
                ax_sns = sns.regplot(
                    data=temp_df_diag, x=resp_name, y=cyto, ax=ax, scatter=False,
                    color=np.unique(ListedColormap(diseasecomb_colors.colors[np.unique(cyto_colorseq_diag)]).colors),
                    label=None, robust=True)
                plt.setp(ax_sns.lines, zorder=1)
                # Add text: Correlation value and p-value per diagnosis
                # x: temp_df.max()[1] / 2 - temp_df.max()[1] / 10, y : temp_df.max()[0] + 1
                ax.text(text_pos[diag], temp_df.max()[0] + 1,
                        r'{}: r = {:.2f}; p = {:.2e}'.format(diag, corr_pval[0][2], corr_pval[0][3]),
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
            save_folder, "{}_diagnosis_{}_transcripts_r{}_{}_{}.{}".format(
                'FigG-I', corr_method, str(distance), cyto, resp_name, fileformat)))
        plt.close()

    return dict_corr
