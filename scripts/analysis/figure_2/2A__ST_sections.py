from scripts.utils import gene_lists, add_observables as ctools

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
        color_dict[cyto] = "#ff7f00"  # orange LICHEN "#ff7f00"
    elif cyto == "IL13":
        color_dict[cyto] = "#e41a1c"  # red AE
    else:
        color_dict[cyto] = "#377eb8"  # blue PSO

    diseases = ["PSO", "AE", "LICHEN", "PRP"]
    biopsy_type = ["LESONAL", "NON LESIONAL"]

    # samples = np.unique(sup_adata[cyto].obs['sample'])
    samples, crops_img = get_cropped_sampleimg(img_key=img_key)
    for ind, sample in enumerate(samples):
        temp_adata = sup_adata[sup_adata.obs["sample"] == sample]
        # read out counts of cytokine and responder spots
        test_counts = temp_adata[temp_adata.obs[obs_label] == cyto]
        # get labels of cytokine and responder spots
        cell_types_unique = list(np.unique(temp_adata.obs[obs_label]))

        if len(cell_types_unique) > 0:

            diagnose = temp_adata.obs[diseases].loc[:, (temp_adata.obs[diseases] != 0).all()].columns.values[0]
            fig, ax = plt.subplots(figsize=fig_size)
            sc.pl.spatial(test_counts, size=1.1, img_key=img_key, library_id=sample,
                          color=obs_label, alpha_img=1, show=False, crop_coord=crops_img[ind],
                          ax=ax, palette=['black'], legend_loc='no')
            sc.pl.spatial(test_counts, size=size, img_key=img_key, library_id=sample,
                          color=obs_counts, alpha_img=0, show=False, crop_coord=crops_img[ind],
                          ax=ax, vmin=1,
                          title=" ".join(["Diagnose:", diagnose, "; Biopsy type:",
                                          temp_adata.obs[biopsy_type].loc[:, (temp_adata.obs[biopsy_type] !=
                                                                              0).all()].columns.values[0]]))
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

        # plot_cytokinecounts(adata, cyto, save_folder, obs_label=obs_label, obs_counts=obs_counts, img_key=img_key)

        # Get max. UMI-counts per tissue section
        samples = np.unique(adata.obs['sample'].values)
        max_umicounts = []
        for sample in samples:
            sample_adata = adata[adata.obs['sample'] == sample].copy()
            m_cyto_counts = sample_adata.obs["_".join([cyto, 'clusters'])] == cyto
            max_umicounts.append(sample_adata.obs["_".join([cyto, 'counts'])][m_cyto_counts].sum())

        print("Max. UMI-counts in tissue section for {}: ".format(cyto), np.amax(max_umicounts))
        print("Mean. UMI-counts in tissue section for {}: ".format(cyto), np.mean(max_umicounts))

        if cyto == 'IL17A':
            diagnose = 'PSO'
        elif cyto == 'IL13':
            diagnose = 'AE'
        elif cyto == 'IFNG':
            diagnose = 'LICHEN'
        umicounts = []
        for sample in samples:
            sample_adata = adata[np.all([adata.obs['sample'] == sample,
                                         adata.obs['disease'] == diagnose], axis=0)].copy()
            m_cyto_counts = sample_adata.obs["_".join([cyto, 'clusters'])] == cyto
            umicounts.append(sample_adata.obs["_".join([cyto, 'counts'])][m_cyto_counts].sum())

        print("Max. UMI-counts in disease for {}: ".format(cyto), np.amax(umicounts))
        print("Mean. UMI-counts in disease for {}: ".format(cyto), np.mean(umicounts))

    return adata


def main(save_folder, adata, img_key):
    # 1. Get cytokines and responders
    t_cell_cytocines, cyto_resps_list, cytokine_responders = gene_lists.get_publication_cyto_resps()

    # 3. Add meta data like which samples belong to which donor (optional)
    if "patient" not in adata.obs_keys():
        adata, tissue_cell_labels, disease_labels, lesion_labels = ctools.add_metadata(adata)
        # 1.2 Remove spots having no tissue/cell labels (since 06.10.2020)
        adata = adata[np.where(adata.obs[tissue_cell_labels].to_numpy().any(axis=1))[0]]

    """Paper Figure 4B: Highlight cytokine and responder genes containing spots and UMI-counts """
    convert_categories_cytokines_responders_others(adata, cyto_responder_genes=cytokine_responders,
                                                   save_folder=save_folder, img_key=img_key)


if __name__ == '__main__':
    today = date.today()
    # create saving folder in current project path
    savepath = os.path.join(os.environ['PYTHONPATH'].split(os.pathsep)[0], "output", "Figure_2A", str(today))
    os.makedirs(savepath, exist_ok=True)

    # 2. Load unpre-processed anndata object
    unpp_st_adata = sc.read(os.path.join(os.environ['PYTHONPATH'].split(os.pathsep)[0],
                                         "adata_storage/2020-10-06/st_adata_P15509_P16357_wo_4_7_unpp.h5"))

    img_resolution = 'hires'
    main(save_folder=savepath, adata=unpp_st_adata, img_key=img_resolution)
