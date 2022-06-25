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
# sc.set_figure_params(color_map='viridis')
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

    samples, crops_img = helper_functions.get_cropped_sampleimg(img_key=img_key)
    # for ind, specimen in enumerate(sup_adata.obs["specimen"].cat.categories):
    for ind, specimen in enumerate(samples):
        # temp_adata = sup_adata[sup_adata.obs["specimen"] == specimen]
        temp_adata = sup_adata[sup_adata.obs["sample"] == specimen]
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

            fig, ax = plt.subplots(figsize=fig_size, ncols=1, constrained_layout=True,
                                   gridspec_kw={'width_ratios': [10]})
            # crop_coord=crops_img[ind]
            sc.pl.spatial(temp_adata, size=1.3, img_key=img_key, library_id=library_id_temp,
                          color=obs_label, groups=cell_types_unique, zorder=1, show=False,
                          legend_loc='left margin', ax=ax, alpha_img=1, crop_coord=crops_img[ind],
                          palette=list_colors, title='')
            if cyto == "IL17A":
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=library_id_temp,
                              color=obs_counts, zorder=2, alpha_img=0, show=False, title='',
                              ax=ax, vmin=0, vmax=800, crop_coord=crops_img[ind])
            elif cyto == "IFNG":
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=library_id_temp,
                              color=obs_counts, zorder=2, alpha_img=0, show=False, title='',
                              ax=ax, vmin=0, vmax=100, crop_coord=crops_img[ind])
            elif cyto == "IL13":
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=library_id_temp,
                              color=obs_counts, zorder=2, alpha_img=0, show=False, title='',
                              ax=ax, vmin=0, vmax=50, crop_coord=crops_img[ind])
            else:
                sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=library_id_temp,
                              color=obs_counts, zorder=2, alpha_img=0, show=False, title='',
                              ax=ax, vmin=0, vmax=100, crop_coord=crops_img[ind])
            # Invert both axis due to flipped and mirrored images
            ax.invert_xaxis()
            ax.invert_yaxis()

            # Colorbar label
            cbar_ax = fig.axes[-1]
            cbar_ax.get_yaxis().labelpad = 15
            cbar_ax.get_yaxis().spacing = 'proportional'
            cbar_ax.get_yaxis().format = '%1i'
            cbar_ax.set_ylabel('UMI-counts / spot', rotation=90, fontsize=legend_fontsize)
            cbar_ax.tick_params(labelsize=xy_ticks)

            # Legend: outside of axis
            leg = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3, fontsize=legend_fontsize)
            leg.get_frame().set_linewidth(0.0)

            # if img_key == "lowres":
            # TODO .pdf still blurry, temp fix with dpi=500
            fig.savefig(os.path.join(
                save_folder, "Radial_plot_{}_{}_{}_{}.pdf".format(specimen, diagnose, biopsy_type_temp, cyto)),
                bbox_inches='tight',  bbox_extra_artists=(leg,), format="pdf", dpi=500)
            plt.close(fig=fig)
            # else:
            #     plt.savefig(os.path.join(save_folder, "_".join(["Radial_plot", specimen, cyto, fileformat])),
            #                 bbox_inches='tight',  bbox_extra_artists=(leg,))
            #     plt.close()
