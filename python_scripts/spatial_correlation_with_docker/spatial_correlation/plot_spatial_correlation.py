"""Visualise Weighted Correlation
    File name: plot_spatial_correlation.py
    Author: Christina Hillig
    Date created: February/01/2021
    Date last modified: September/01/2025
    Python Version: 3.8
"""

# import scripts
import corr_statistics as corr_stats
import utils_plots
import plot_colorbar_legend

# Plotting packages
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

# System specific
import os

# Calculation packages
import scanpy as sc
import numpy as np
import pandas as pd


# Global constants (adjust if needed)
fig_size, title_fontsize, axis_label_fontsize, legend_fontsize, fileformat, img_key, xy_ticks, text_fontsize = \
    utils_plots.figure_params()


def get_linefit(df, resp_name, cyto):
    # Save df in a new df and sort it - otherwise fill_between won't work
    xy_vals = df[[resp_name, cyto]].copy()
    xy_vals = xy_vals.sort_values(by=resp_name)

    # Prepare data for linear regression fit
    x = np.asarray([xy_vals[resp_name].values]).T
    y = np.asarray([xy_vals[cyto].values]).T

    # Fit Line going through origin (0, 0), intercept of 0
    est_fitted, ypred, df_results = corr_stats.create_ls_model(
        x=x, y=y, w=None, const=False)

    return est_fitted, ypred, df_results, x


def plot__stwc_generic(
    df_counts, cytokine_responders, save_folder, distance,
    corr_method, dict_corr, group_col=None, group_colors=None,
    colorbar_fn=None, filename_suffix="",
    color_dots=True, size_col=None, size_multiplier=36, size_legend_fn=None
):
    """
    Generic plotting function for cytokine/responder correlations.

    Parameters
    ----------
    df_counts : pd.DataFrame
    cytokine_responders : list
    save_folder : str
    distance : int
    corr_method : str
    dict_corr : dict
    group_col : str or None, column for grouping/coloring
    group_colors : dict or None, mapping from group labels -> colors or ints
    colorbar_fn : callable or None, function to plot standalone colorbar
    filename_suffix : str, added to output filename
    color_dots : bool, whether to color dots (else black)
    size_col: None,
    size_multiplier: int, =30,
    size_legend_fn: callable or None, function to plot standalone legend of point sizes
    """
    for cyto in cytokine_responders:
        resp_name = f"{cyto}_responder"
        temp_df = df_counts[[cyto, resp_name, "Cluster_size", "Cluster_num_spots"]].copy()

        # Cast numeric
        temp_df = temp_df.dropna(subset=[cyto, resp_name])
        temp_df[cyto] = temp_df[cyto].astype(float)
        temp_df[resp_name] = temp_df[resp_name].astype(float)

        # Fit regression line
        est_fit, ypred, df_fit, x = get_linefit(df=temp_df, resp_name=resp_name, cyto=cyto)

        # Correlation stats
        corr_pval = [a for a in dict_corr[cyto][corr_method] if a[0] == cyto][0]

        # Plot
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)

        # Regression line + error band
        ax.plot(x, ypred, color='black', lw=2, zorder=1)
        ax.fill_between(x.transpose()[0],
                        df_fit['mean_ci_lower'].values,
                        df_fit['mean_ci_upper'].values,
                        alpha=.1, color='black', lw=0.1)

        # Assign colors
        if color_dots and group_col and group_colors:
            cyto_colorseq = [group_colors[val] for val in df_counts.loc[temp_df.index, group_col]]
        else:
            cyto_colorseq = 'black'

        # Assign sizes
        if size_col and "Cluster_size" in temp_df:
            temp_df["Cluster size"] = (temp_df["Cluster_size"].values**2 * size_multiplier).astype(int)
            dot_sizes = temp_df["Cluster size"]
        else:
            dot_sizes = 100  # fixed size

        # Scatter
        ax.scatter(temp_df[resp_name], temp_df[cyto],
                   c=cyto_colorseq, s=dot_sizes,
                   alpha=0.8, edgecolors='k', lw=0.5, zorder=2)

        # Axes
        ax.set_xlabel("cluster proximity gene(s) counts", fontsize=axis_label_fontsize)
        ax.set_ylabel(f"cluster {cyto.split('_')[0]} counts", fontsize=axis_label_fontsize)

        # Add correlation text
        ax.text(1, temp_df[cyto].max() + 1,
                f"weighted {corr_method.capitalize()} correlation = {corr_pval[2]:.2f}\np-value = {corr_pval[3]:.2e}",
                fontstyle='italic', fontsize=10)

        ax.tick_params(labelsize=16)
        sns.despine(ax=ax)
        plt.tight_layout()

        # Save
        fname = os.path.join(
            save_folder,
            f"{filename_suffix}_ColorDots_{color_dots}_Clustersize){size_col}_{corr_method}_r{distance}_{cyto}_{resp_name}{fileformat}"
        )
        fig.savefig(fname, bbox_inches="tight")
        plt.close(fig)

        # === Standalone colorbar legend ===
        if color_dots and group_col and group_colors and colorbar_fn:
            cmap = mpl.colors.ListedColormap(list(group_colors.values()))
            labels = list(group_colors.keys())
            colorbar_fn(cmap, labels, save_folder, group_col)

        # === Standalone size legend ===
        if size_col and size_legend_fn:
            uniq_vals = np.unique(df_counts[size_col].astype(int))
            size_legend_fn(uniq_vals, size_multiplier, save_folder)


# -----------------------------
# Wrappers
# -----------------------------

def plot__stwc_biopsytype(df_counts, cytokine_responders, save_folder,
                          distance, corr_method, dict_corr,
                          color_dots=True, scale_dot_size=True):
    unique_types = df_counts["biopsy_type"].unique()
    colors = dict(zip(unique_types, sc.pl.palettes.default_102[:len(unique_types)]))

    plot__stwc_generic(
        df_counts, cytokine_responders, save_folder, distance, corr_method, dict_corr,
        group_col="biopsy_type", group_colors=colors,
        colorbar_fn=plot_colorbar_legend.plot_standalone_colorbar_obs,
        filename_suffix="biopsytype",
        color_dots=color_dots,
        size_col="Cluster_size" if scale_dot_size else None,
        size_multiplier=36,
        size_legend_fn=plot_colorbar_legend.plot_standalone_legend
    )


def plot__stwc_disease(df_counts, cytokine_responders, save_folder,
                       distance, corr_method, dict_corr,
                       color_dots=True, scale_dot_size=True):
    unique_disease = df_counts["disease"].unique()
    colors = dict(zip(unique_disease, sc.pl.palettes.default_102[:len(unique_disease)]))

    plot__stwc_generic(
        df_counts, cytokine_responders, save_folder, distance, corr_method, dict_corr,
        group_col="disease", group_colors=colors,
        colorbar_fn=plot_colorbar_legend.plot_standalone_colorbar_obs,
        filename_suffix="disease",
        color_dots=color_dots,
        size_col="Cluster_size" if scale_dot_size else None,
        size_multiplier=36,
        size_legend_fn=plot_colorbar_legend.plot_standalone_legend
    )


def plot__stwc_tissuelayers(df_counts, cytokine_responders, save_folder,
                            distance, corr_method, dict_corr,
                            color_dots=True, scale_dot_size=True):
    unique_layers = df_counts["tissue_layer"].unique()
    colors = dict(zip(unique_layers, sc.pl.palettes.default_102[:len(unique_layers)]))

    plot__stwc_generic(
        df_counts, cytokine_responders, save_folder, distance, corr_method, dict_corr,
        group_col="tissue_layer", group_colors=colors,
        colorbar_fn=plot_colorbar_legend.plot_standalone_colorbar_obs,
        filename_suffix="tissuelayers",
        color_dots=color_dots,
        size_col="Cluster_size" if scale_dot_size else None,
        size_multiplier=36,
        size_legend_fn=plot_colorbar_legend.plot_standalone_legend
    )


def plot__stwc_patients(df_counts, cytokine_responders, save_folder,
                        distance, corr_method, dict_corr,
                        color_dots=True, scale_dot_size=True):
    unique_patients = df_counts["Patient"].astype(int).unique()
    colors = dict(zip(unique_patients, sc.pl.palettes.default_102[:len(unique_patients)]))

    plot__stwc_generic(
        df_counts, cytokine_responders, save_folder, distance, corr_method, dict_corr,
        group_col="Patient", group_colors=colors,
        colorbar_fn=plot_colorbar_legend.plot_standalone_colorbar_obs,
        filename_suffix="patients",
        color_dots=color_dots,
        size_col="Cluster_size" if scale_dot_size else None,
        size_multiplier=36,
        size_legend_fn=plot_colorbar_legend.plot_standalone_legend
    )
