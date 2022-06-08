#! /usr/bin/python
"""Plot counts of cytokine positive spots on H&E slides
    File name: Fig2A__ST_scetions.py
    Author: Christina Hillig
    Date created: October/xx/2020
    Date last modified: May/02/2021
    Python Version: 3.7
"""
from python_scripts.utils import gene_lists

import scanpy as sc
import numpy as np
import os
from datetime import date


from collections import OrderedDict
import matplotlib.pyplot as plt

# figure properties
fig_size = (8, 8)
xy_fontsize = 16
xy_ticks = 12
title_fontsize = 18
legend_fontsize = 14
text_fontsize = 14
fileformat = '.pdf'


def get_cropped_sampleimg(img_key):
    samples = ['P15509_1003', 'P15509_1004', 'P16357_1003', 'P16357_1019', 'P16357_1020', 'P16357_1031',
               'P16357_1032', 'P16357_1036']
    if img_key == 'lowres':
        crops_img = [(100, 400, 520, 130), (49, 349, 510, 100), (80, 380, 510, 100), (150, 450, 510, 100),
                     (150, 450, 530, 120), (135, 435, 530, 120), (150, 450, 510, 100), (160, 460, 510, 100)]
    else:
        res_proportion = 2000 / 600
        crops_img = [(100, 400, 520, 130), (49, 349, 510, 100), (80, 380, 510, 100), (150, 450, 510, 100),
                     (150, 450, 530, 120), (135, 435, 530, 120), (150, 450, 510, 100), (160, 460, 510, 100)]
        crops_img = (np.asarray(crops_img) * res_proportion).round()

    return samples, crops_img


def plot_cytokinecounts(sup_adata, cyto, save_folder, obs_label, obs_counts, img_key="lowres"):
    """Plot counts of cytokines on H&E image

    Parameters
    ----------
    sup_adata : annData
        containing annData object
    cyto : str
        cytokine name
    save_folder : str
    obs_label : str
    obs_counts : str
    img_key : str
        hires or lowres

    Returns
    -------

    """
    # Size of count spots
    size = 0.9

    color_dict = OrderedDict()
    if cyto == 'IFNG':
        color_dict[cyto] = "#ff7f00"  # orange LP
    elif cyto == "IL13":
        color_dict[cyto] = "#e41a1c"  # red AD
    else:
        color_dict[cyto] = "#377eb8"  # blue Pso

    # samples = np.unique(sup_adata[cyto].obs['sample'])
    samples, crops_img = get_cropped_sampleimg(img_key=img_key)
    for ind, sample in enumerate(samples):
        temp_adata = sup_adata[sup_adata.obs["sample"] == sample].copy()
        # read out counts of cytokine and responder spots
        test_counts = temp_adata[temp_adata.obs[obs_label] == cyto].copy()
        # get labels of cytokine and responder spots
        # cell_types_unique = list(np.unique(temp_adata.obs[obs_label]))

        if test_counts.shape[0] > 0:
            diagnose = " & ".join(test_counts.obs['DISEASE'].cat.categories)
            biopsy_type_temp = " & ".join(test_counts.obs['biopsy_type'].cat.categories)
            library_id_temp = "_".join(test_counts.obs['library_id'].cat.categories[0].split('_')[:-1])

            fig, ax = plt.subplots(figsize=fig_size)
            sc.pl.spatial(test_counts, size=1.1, img_key=img_key, library_id=library_id_temp,
                          color=obs_label, alpha_img=1, show=False, crop_coord=crops_img[ind],
                          ax=ax, palette=['black'], legend_loc='no')
            sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=library_id_temp,
                          color=obs_counts, alpha_img=0, show=False, crop_coord=crops_img[ind],
                          ax=ax, vmin=1,
                          vmax=sup_adata[sup_adata.obs[obs_label] == cyto].obs["{}_counts".format(cyto)].max(),
                          title=" ".join(["Diagnose:", diagnose, "; Biopsy type:", biopsy_type_temp]))
            # Invert both axis due to flipped and mirrored images
            ax.invert_xaxis()
            ax.invert_yaxis()

            # Colorbar label
            cbar_ax = fig.axes[-1]
            cbar_ax.get_yaxis().labelpad = 12
            cbar_ax.set_ylabel('UMI-counts', rotation=90, fontsize=legend_fontsize)
            cbar_ax.tick_params(labelsize=xy_ticks)

            plt.tight_layout()

            if img_key == "lowres":
                plt.savefig(os.path.join(save_folder, "_".join(["Radial_plot", sample, cyto, ".png"])),
                            bbox_inches='tight')
                plt.close()
            else:
                plt.savefig(os.path.join(save_folder, "_".join(["Radial_plot", sample, cyto, fileformat])),
                            bbox_inches='tight')
                plt.close()


def convert_categories_cytokines_responders_others(adata, cyto_responder_genes, save_folder, img_key):
    """Add observable to annData object containing the three categories cytokine, responders and others

    Parameters
    ----------
    adata : annData
    cyto_responder_genes : dict
    save_folder : str
    img_key : str

    Returns
    -------

    """

    for cyto in cyto_responder_genes.keys():
        obs_counts = "_".join([cyto, 'counts'])
        obs_label = "_".join([cyto, 'responders'])
        obs_clusters = "_".join([cyto, 'clusters'])

        # create array with others
        if "counts" in adata.layers.keys():
            default_counts = np.copy(adata.layers['counts']).sum(axis=1)
        else:
            default_counts = np.copy(adata.X).sum(axis=1)
        default_label = np.array(['Others'] * adata.shape[0], dtype='<U10')
        default_clusters = np.array([2] * adata.shape[0])

        # Get available responder genes
        available_responders = list(set(adata.var.index) & set(cyto_responder_genes[cyto]))

        # 2. Get mask for responder genes spot positions
        # 2.1 Find index of responder genes
        varindex_responder_genes = np.where(
            adata.var.index[np.newaxis, :] == np.array(available_responders)[:, np.newaxis])[1]
        # 2.2 Get count matrix of responder genes
        if "counts" in adata.layers.keys():
            responder_matrix = np.copy(adata.layers['counts'])[:, varindex_responder_genes]
        else:
            responder_matrix = np.copy(adata.X)[:, varindex_responder_genes]
        # Identify where counts > 0
        m_responders = responder_matrix > 0
        m_array_responders = m_responders.sum(axis=1).astype(bool)

        # Add label of responder
        default_label[m_array_responders] = 'Responders'
        default_clusters[m_array_responders] = 1
        default_counts[m_array_responders] = np.sum(responder_matrix, axis=1)[m_array_responders]

        # 3. Get cytokine spot positions
        if cyto in adata.var_names:
            if "counts" in adata.layers.keys():
                cyto_matrix = np.copy(adata.layers['counts'])[:, np.where(adata.var.index == cyto)[0]]
            else:
                cyto_matrix = np.copy(adata.X)[:, np.where(adata.var.index == cyto)[0]]
            m_cyto = cyto_matrix[:, 0] > 0

            # Add label cytokine positive spots
            default_clusters[m_cyto] = 0
            default_label[m_cyto] = cyto
            default_counts[m_cyto] = cyto_matrix[:, 0][m_cyto]

        # 4. Add observables to adata
        adata.obs[obs_clusters] = default_clusters
        adata.obs[obs_label] = default_label
        adata.obs[obs_counts] = default_counts

        # 5. convert to categorical
        adata.obs[obs_clusters] = adata.obs[obs_label].astype('category')
        adata.obs[obs_label] = adata.obs[obs_label].astype('category')

        # Plot cytokine counts on H&E image
        plot_cytokinecounts(adata, cyto, save_folder, obs_label=obs_label, obs_counts=obs_counts, img_key=img_key)

        # Get max. UMI-counts per tissue section
        specimens = np.unique(adata.obs['specimen'].values)
        max_umicounts = []
        for specimen in specimens:
            sample_adata = adata[adata.obs['specimen'] == specimen].copy()
            m_cyto_counts = sample_adata.obs["_".join([cyto, 'clusters'])] == cyto
            max_umicounts.append(sample_adata.obs["_".join([cyto, 'counts'])][m_cyto_counts].sum())

        print("Max. UMI-counts in tissue section for {}: ".format(cyto), np.amax(max_umicounts))
        print("Mean. UMI-counts in tissue section for {}: ".format(cyto), np.mean(max_umicounts))

        if cyto == 'IL17A':
            diagnose = 'Pso'
        elif cyto == 'IL13':
            diagnose = 'AD'
        elif cyto == 'IFNG':
            diagnose = 'LP'
        else:
            diagnose = 'Unknown'
        umicounts = []
        for specimen in specimens:
            sample_adata = adata[np.all([adata.obs['specimen'] == specimen,
                                         adata.obs['DISEASE'] == diagnose], axis=0)].copy()
            m_cyto_counts = sample_adata.obs["_".join([cyto, 'clusters'])] == cyto
            umicounts.append(sample_adata.obs["_".join([cyto, 'counts'])][m_cyto_counts].sum())

        print("Max. UMI-counts in disease for {}: ".format(diagnose), np.amax(umicounts))
        print("Mean. UMI-counts in disease for {}: ".format(diagnose), np.mean(umicounts))

    return adata


def main(save_folder, adata):
    img_key = 'hires'
    # 1. Get cytokines and responders
    t_cell_cytocines, cyto_resps_list, cytokine_responders = gene_lists.get_publication_cyto_resps()

    """Paper Figure 4B: Highlight cytokine and responder genes containing spots and UMI-counts """
    convert_categories_cytokines_responders_others(adata, cyto_responder_genes=cytokine_responders,
                                                   save_folder=save_folder, img_key=img_key)


if __name__ == '__main__':
    today = date.today()
    # create saving folder in current project path
    savepath = os.path.join("..", "..", "..", "output", "Figure_2A", str(today))
    os.makedirs(savepath, exist_ok=True)

    # 2. Load unpre-processed anndata object
    unpp_st_adata = sc.read(
        os.path.join("..", "..", "..", "adata_storage", "2022-04-08",
                     "Spatial Transcriptomics_unpp_cleaned_PsoADLP.h5"))

    img_resolution = 'hires'
    main(save_folder=savepath, adata=unpp_st_adata)
