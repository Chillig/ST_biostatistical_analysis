"""Visualise count distributions of responder and cytokine genes
    File name: plot_count_distributions.py
    Author: Christina Hillig
    Date created: May/01/2021
    Date last modified: May/01/2021
    Python Version: 3.7
"""

# import scripts
from python_scripts.spatial_correlation import corr_statistics as corr_stats
from python_scripts.spatial_correlation import helper_functions

# Plotting packages
import matplotlib.pyplot as plt
import seaborn as sns

# System specific
import os

# Calculation packages
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.stats as scstats


# Figure params
sc.set_figure_params(color_map='viridis')
fig_size, title_fontsize, axis_label_fontsize, legend_fontsize, fileformat, img_key, xy_ticks, text_fontsize = \
    helper_functions.figure_params()


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
