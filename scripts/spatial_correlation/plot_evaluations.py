"""Evaluation metrics
    File name: plot_evaluations.py
    Author: Christina Hillig
    Date created: May/01/2021
    Date last modified: May/01/2021
    Python Version: 3.7
"""

# import scripts
from scripts.spatial_correlation import helper_functions

# Plotting packages
import matplotlib.pyplot as plt

# System specific
import os

# Calculation packages
import scanpy as sc
import numpy as np


# Figure params
sc.set_figure_params(color_map='viridis')
fig_size, title_fontsize, axis_label_fontsize, legend_fontsize, fileformat, img_key, xy_ticks, text_fontsize = \
    helper_functions.figure_params()


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
