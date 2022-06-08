"""Plot conditional density clusters and graphs with and without H&E image
    File name: plot_clusters.py
    Author: Christina Hillig
    Date created: May/01/2021
    Date last modified: May/01/2021
    Python Version: 3.7
"""

# import scripts
from python_scripts.spatial_correlation import helper_functions

# Plotting packages
import matplotlib.pyplot as plt

# System specific
import os

# Calculation packages
from collections import OrderedDict
import scanpy as sc
import numpy as np


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
        color_dict[cyto] = "#ff7f00"  # orange LICHEN
    elif cyto == "IL13":
        color_dict[cyto] = "#e41a1c"  # red AE
    else:
        color_dict[cyto] = "#377eb8"  # blue PSO

    color_dict["Responders"] = "y"
    color_dict["Others"] = "silver"

    # samples, crops_img = helper_functions.get_cropped_sampleimg(img_key=img_key)
    for ind, specimen in enumerate(sup_adata.obs["specimen"].cat.categories):
        temp_adata = sup_adata[sup_adata.obs["specimen"] == specimen]
        # read out counts of cytokine and responder spots
        test_counts = temp_adata[temp_adata.obs[obs_label] != "Others"]
        # get labels of cytokine and responder spots
        cell_types_unique = list(np.unique(temp_adata.obs[obs_label]))

        list_colors = []
        for i_color in cell_types_unique:
            list_colors.append(color_dict[i_color])

        if test_counts.shape[0] > 0:
            diagnose = " & ".join(test_counts.obs['DISEASE'].cat.categories)
            biopsy_type_temp = " & ".join(test_counts.obs['biopsy_type'].cat.categories)
            library_id_temp = "_".join(test_counts.obs['library_id'].cat.categories[0].split('_')[:-1])

            fig, ax = plt.subplots(figsize=fig_size)
            # crop_coord=crops_img[ind]
            sc.pl.spatial(temp_adata, size=1.3, img_key=img_key, library_id=library_id_temp,
                          color=obs_label, groups=cell_types_unique, zorder=1, show=False,
                          legend_loc='left margin', ax=ax, alpha_img=1,
                          palette=list_colors)
            if cyto == "IL17A":
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=library_id_temp,
                              color=obs_counts, zorder=2, alpha_img=0, show=False,
                              ax=ax, vmin=0, vmax=800,
                              title=" ".join(["Diagnose:", diagnose, "; Biopsy type:", biopsy_type_temp]))
            elif cyto == "IFNG":
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=library_id_temp,
                              color=obs_counts, zorder=2, alpha_img=0, show=False,
                              ax=ax, vmin=0, vmax=100,
                              title=" ".join(["Diagnose:", diagnose, "; Biopsy type:", biopsy_type_temp]))
            elif cyto == "IL13":
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=library_id_temp,
                              color=obs_counts, zorder=2, alpha_img=0, show=False,
                              ax=ax, vmin=0, vmax=50,
                              title=" ".join(["Diagnose:", diagnose, "; Biopsy type:", biopsy_type_temp]))
            else:
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=library_id_temp,
                              color=obs_counts, zorder=2, alpha_img=0, show=False,
                              ax=ax, vmin=0, vmax=100,
                              title=" ".join(["Diagnose:", diagnose, "; Biopsy type:", biopsy_type_temp]))
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
                plt.savefig(os.path.join(save_folder, "_".join(["Radial_plot", specimen, cyto, ".png"])),
                            bbox_inches='tight',  bbox_extra_artists=(leg,))
                plt.close()
            else:
                plt.savefig(os.path.join(save_folder, "_".join(["Radial_plot", specimen, cyto, fileformat])),
                            bbox_inches='tight',  bbox_extra_artists=(leg,))
                plt.close()