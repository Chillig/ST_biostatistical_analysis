# import scripts
from python_scripts.spatial_correlation import corr_statistics as corr_stats
from python_scripts.spatial_correlation import helper_functions
from python_scripts.spatial_correlation.plots import plot_colorbar_legend

import pandas as pd
import numpy as np
import os

# Plotting packages
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns

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


def plot__stwc_tissuelayers(bulk_df_counts, st_df_counts, cytokine_responders, save_folder, distance, corr_method):
    """Calculate Correlation between each Cytokine positive spot and its nn responder genes spots for each sample
    -> Plot for Figure 4G-I

    Parameters
    ----------
    bulk_df_counts : pandas.Dataframe
    st_df_counts : pandas.Dataframe
    cytokine_responders : dict
    save_folder : str
    distance : int
    corr_method : str

    Returns
    -------

    """
    # Get unique combinations of the tissue layer to assign to each a color before the plotting
    # 1. map each string to an integer value
    color_seq = pd.factorize(st_df_counts['tissue_layer'])[0]
    num_unique_colors = len(np.unique(color_seq))
    # 2. assign a color to each combination
    dict_tissue_color = OrderedDict(zip(st_df_counts['tissue_layer'], color_seq))
    tissuecomb_colors = plt.cm.get_cmap("tab20", num_unique_colors)

    # 3. Plot Colorbar with tissue layer combinations
    plot_colorbar_legend.plot_standalone_colorbar(
        tissuecomb_colors=tissuecomb_colors, labels=dict_tissue_color.keys(), save_folder=save_folder)

    # 4. get unique cluster sizes
    cluster_sizes = np.unique(st_df_counts['Cluster_size'])
    plot_colorbar_legend.plot_standalone_legend(cluster_sizes, save_folder=save_folder, size_multiplier=size_multiplier)

    # 5. BULK: convert nan values to 0
    bulk_df_counts = bulk_df_counts.replace({np.nan: 0})

    dict_corr_st = dict({'pearson': [], 'spearman': []})
    dict_corr_bulk = dict({'pearson': [], 'spearman': []})
    for cyto in cytokine_responders:
        resp_name = "_".join([cyto, 'responder'])
        temp_df = st_df_counts[[cyto, resp_name, "Cluster_size", "disease", 'tissue_layer']].copy()
        weighted_cytoname = "_".join(['weighted', cyto])

        temp_df = temp_df[~np.isnan(temp_df[cyto].values.astype(np.float64))]
        temp_df[cyto] = temp_df[cyto].values.astype(np.float)
        temp_df[resp_name] = temp_df[resp_name].values.astype(np.float)

        # Read out Pseudo-bulk approach data
        bulk_df_temp = bulk_df_counts[[cyto, resp_name, 'disease', weighted_cytoname]]

        # # map each string to an integer value
        cyto_colorseq = []
        for ind_dict, value in enumerate(temp_df['tissue_layer']):
            cyto_colorseq.append(dict_tissue_color[value])

        # 1. Calculate correlation
        # 1.1 Apply weights to cytokine counts and add it to df
        weighted_cytocounts = temp_df[cyto]
        temp_df[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)

        dict_corr_st = get_correlation_stats(
            df=temp_df, resp_name=resp_name, cyto=cyto, weight_name=weighted_cytoname, dict_corr=dict_corr_st)
        dict_corr_bulk = get_correlation_stats(
            df=bulk_df_temp, resp_name=resp_name, cyto=cyto, weight_name=weighted_cytoname, dict_corr=dict_corr_bulk)

        get_sizes = np.unique(temp_df['Cluster_size'])
        if any(get_sizes > 6):
            increased_ylimit = 3
        else:
            increased_ylimit = 1

        # Increase point size
        temp_df["Cluster size"] = (temp_df["Cluster_size"].values**2 * size_multiplier).astype(int)

        # Fit Line going through origin (0, 0), intercept of 0
        est_woweights_fitted, ypred_woweights, df_woweights, x_st = get_linefit(
            df=temp_df, resp_name=resp_name, cyto=cyto)
        bulk_est_woweights_fitted, ypred_woweights_bulk, df_woweights_bulk, x_bulk = get_linefit(
            df=bulk_df_temp, resp_name=resp_name, cyto=cyto)

        #  2. Plot correlation
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        # Unweighted: Plot weighted line fit
        ax.plot(x_st, ypred_woweights, label="ST fit", color='black', zorder=1, lw=2)
        # - ST: plotting error band
        ax.fill_between(x_st.transpose()[0], df_woweights['mean_ci_lower'].values, df_woweights['mean_ci_upper'].values,
                        alpha=.1, label='5 - sigma interval', color='black', lw=0.1)
        # Unweighted: Plot weighted line fit of bulk RNA-seq data, add dashed fitted line
        ax.plot(x_bulk, ypred_woweights_bulk, label="Bulk fit", color='black', zorder=1, lw=2, ls='--')
        # - Bulk: plotting error band
        ax.fill_between(
            x_bulk.transpose()[0], df_woweights_bulk['mean_ci_lower'].values, df_woweights_bulk['mean_ci_upper'].values,
            alpha=.1, label='5 - sigma interval', color='black', lw=0.1)
        # Unweighted: Plot points
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
        corr_pval = [a_tuple for a_tuple in dict_corr_st[corr_method] if a_tuple[0] == cyto]

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

        plt.tight_layout()

        # remove upper and right edge lines in plot
        sns.despine(ax=ax)

        # 3. Save figure -> creates figure for Workflow figure 1
        fig.savefig(os.path.join(
            save_folder, "{}_tissuelayers_weighted_transcripts_{}_r{}_{}_{}{}".format(
                'FigF-H', corr_method, str(distance), cyto, resp_name, fileformat)),)
        plt.close()
