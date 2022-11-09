"""Plot Colorbar and Legend separately
    File name: plot_colorbar_legend.py
    Author: Christina Hillig
    Date created: May/01/2021
    Date last modified: May/01/2021
    Python Version: 3.7
"""

# import scripts
from python_scripts.spatial_correlation import helper_functions

# Plotting packages
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# System specific
import os

# Calculation packages
import scanpy as sc
import numpy as np


# Figure params
sc.set_figure_params(color_map='viridis')
fig_size, title_fontsize, axis_label_fontsize, legend_fontsize, fileformat, img_key, xy_ticks, text_fontsize = \
    helper_functions.figure_params()


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

    fig = plt.figure(figsize=(8, 6))
    # [left, bottom, width, height]
    ax = fig.add_axes([0.1, 0.05, 0.16, 0.9])  # work only if rotation = 0
    # PyCharm bug: https://stackoverflow.com/questions/23248017/cannot-find-reference-xxx-in-init-py-python-pycharm
    cbar = mpl.colorbar.ColorbarBase(
        ax, cmap=ListedColormap(tissuecomb_colors.colors), norm=norm, orientation='vertical')
    cbar.set_ticks(np.linspace(0, num_colors, num_colors + 1)[:-1] + 0.5)
    cbar.ax.set_yticklabels(labels, fontsize=xy_ticks, rotation=0, ha="left")
    ax.set_title('Tissue layers', fontsize=title_fontsize)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "Tissue_layer_colorbar.pdf"))
    plt.close()


def plot_standalone_colorbar_patient(patient_colors, labels, save_folder):
    """Plot a standalone vertical and horizontal Colorbar

    Parameters
    ----------
    patient_colors : matplotlib.colors.ListedColormap
    labels : list of str
    save_folder : str

    Returns
    -------

    """

    num_colors = len(patient_colors.colors)
    norm = mpl.colors.Normalize(vmin=0, vmax=num_colors)

    fig = plt.figure(figsize=(4, 6))
    # [left, bottom, width, height]
    ax = fig.add_axes([0.1, 0.05, 0.16, 0.9])  # work only if rotation = 0
    # PyCharm bug: https://stackoverflow.com/questions/23248017/cannot-find-reference-xxx-in-init-py-python-pycharm
    cbar = mpl.colorbar.ColorbarBase(
        ax, cmap=ListedColormap(patient_colors.colors), norm=norm, orientation='vertical')
    cbar.set_ticks(np.linspace(0, num_colors, num_colors + 1)[:-1] + 0.5)
    cbar.ax.set_yticklabels(labels, fontsize=xy_ticks, rotation=0, ha="left")
    ax.set_title('Patient', fontsize=title_fontsize)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "Patient_colorbar.pdf"))
    plt.close()


def plot_standalone_legend(points, size_multiplier, save_folder):
    """Plot a legend without a graph

    Parameters
    ----------
    points : numpy.array
    size_multiplier: int
    save_folder : str

    Returns
    -------

    """
    # increase point size
    if 0 in points:
        # have to start at 1, otherwise points disappear in plot due to zero size
        large_points = points.astype(int)**2 * size_multiplier + 1
    else:
        large_points = points.astype(int) ** 2 * size_multiplier
    colors = ['k'] * len(points)

    # Draw circles
    patches = [plt.scatter([], [], marker="o", s=large_points.astype(int)[i], color=colors[i],
                           label="{:s}".format(str(points[i]))) for i in range(len(points))]
    plt.close()

    if any(points >= 20):
        label_space = 7
    elif any((points >= 10) | (points < 20)):
        label_space = 5
    elif any((points >= 6) | (points < 10)):
        label_space = 3
    else:
        label_space = 1

    # Create figure without an axis
    fig = plt.figure(figsize=(8, 10))
    fig.legend(patches, points, labelspacing=label_space, title="# cyto+ spots in cluster",
               loc='center', frameon=False, facecolor="w", handletextpad=2, handlelength=2,
               title_fontsize=title_fontsize, fontsize=legend_fontsize,
               bbox_transform=fig.transFigure, borderpad=.0, scatterpoints=1, ncol=3)
    plt.tight_layout()
    # Save figure
    fig.savefig(os.path.join(save_folder, "_".join(["Legend", fileformat])),)
    plt.close()
