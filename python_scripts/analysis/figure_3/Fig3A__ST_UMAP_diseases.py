#! /usr/bin/python
"""Plot UMAP of Spatial Transcriptomics data together with IL-17A positive spots
    File name: Fig3A__ST_UMAP_IL17A_tissuelayers.py
    Author: Christina Hillig
    Date created: December/xx/2020
    Date last modified: May/02/2021
    Python Version: 3.7
"""
from python_scripts.utils import gene_lists, add_observables, get_condition_spots

import scanpy as sc
import numpy as np
import os
from datetime import date
from collections import OrderedDict
import matplotlib.pyplot as plt

import pandas as pd


signatures = OrderedDict()
# publication
signatures["IFNG"] = "#ff7f00"  # orange LICHEN
signatures["IL13"] = "#e41a1c"  # red AE
signatures["IL17A"] = "#377eb8"  # blue PSO
signatures["GAPDH"] = '#4daf4a'  # green GAPDH
signatures["CD2"] = 'tab:cyan'  # cyan
signatures["Cytokines"] = 'purple'  # purple

# Figure params
figure_size = (8, 8)
xy_fontsize = 16
xy_ticks = 12
title_fontsize = 18
legend_fontsize = 14
text_fontsize = 12
file_format = '.pdf'
img_key = 'hires'


def plot_disease(adata, obs_name, title, save_folder):
    """Plot UMAP embedded tissue regions with IL17A on top

    Parameters
    ----------
    adata : annData
    obs_name : str
    title : str
    save_folder : str

    Returns
    -------

    """

    adata.obs[obs_name] = adata.obs[obs_name].astype('category')
    adata.obs[obs_name] = adata.obs[obs_name].cat.reorder_categories(['LP', 'AD', 'Pso'])
    palette = ['#ff7f00', '#e41a1c', '#377eb8']

    fig, ax = plt.subplots(figsize=figure_size)
    sc.pl.umap(adata, color=obs_name, use_raw=True, ax=ax, wspace=0.4, show=False,
               size=20, frameon=True, facecolor='white', palette=palette,
               title="")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    plt.savefig(os.path.join(save_folder, "_".join(['UMAP', title, "DISEASE", file_format])))
    plt.close()


def get_celltypes_data(adata, genes):
    """
    Get adata object containing only those cells/spots which express genes of interest

    :param adata: [annData]
    :param genes: [list or string]
    :return:
    """

    if isinstance(genes, list):
        varindex_cyto_genes = \
             np.where(adata.var.index[np.newaxis, :] == np.array(genes)[:, np.newaxis])[1]
        counts_cyto = adata.layers["counts"][:, varindex_cyto_genes]
    else:
        varindex_cyto_genes = np.where(adata.var.index == genes)[1]
        counts_cyto = adata.layers["counts"][:, varindex_cyto_genes][:, 0]
    # create mask
    m_cyto = counts_cyto > 0
    m_cyto = np.any(m_cyto, axis=1)
    adata = adata.copy()[m_cyto]

    return adata


def get_tissueregions(adata, tissue_label):
    """Merge Epidermis layers into EPIDERMIS and Dermis layers into DERMIS

    :param adata:
    :param tissue_label:
    :return:
    """
    m_epidermis = np.array(
        adata.obs[tissue_label].values)[:, np.newaxis] == np.array(
        ['upper EPIDERMIS', 'basal EPIDERMIS', 'middle EPIDERMIS', 'JUNCTION'])[np.newaxis, :]
    m_epidermis = m_epidermis.sum(axis=1).astype(bool)

    m_dermis = np.array(
        adata.obs[tissue_label].values)[:, np.newaxis] == np.array(
        ['DERdepth1', 'DERdepth2', 'DERdepth3', 'DERdepth4', 'DERdepth5', 'DERdepth6', 'DERdepth7'])[np.newaxis, :]
    m_dermis = m_dermis.sum(axis=1).astype(bool)

    adata.obs['tissue_regions'] = 'Unknown_tissueregions'
    adata.obs['tissue_regions'][m_dermis] = 'Dermis'
    adata.obs['tissue_regions'][m_epidermis] = 'Epidermis'
    adata.obs['tissue_regions'] = adata.obs['tissue_regions'].astype('category')

    return adata


def main(save_folder, spatial_adata):
    """Read out data for ST DGE Analysis and create UMAPs for Figure 3A

    :return:
    """
    spatial_cluster_label = 'tissue_layer'

    # load data
    cytokines, allinone, cytoresps_dict = gene_lists.get_publication_cyto_resps()
    leukocyte_markers = gene_lists.leukocyte_markers()

    # remove all spots without a tissue label
    spatial_adata = spatial_adata[spatial_adata.obs[spatial_cluster_label] != 'Unknown'].copy()

    # 1. get observable for cytokine genes
    spatial_adata, obs_name = add_observables.convert_variable_to_observable(
        adata=spatial_adata, gene_names=cytokines, task='cell_gene', label='celltype', condition=None)

    spatial_adata, _ = add_observables.convert_variable_to_observable(
        adata=spatial_adata, gene_names=leukocyte_markers,
        task='cell_gene', label='celltype', condition=None)

    # 2. Highlight diseases
    spatial_adata = spatial_adata[spatial_adata.obs['biopsy_type'] == 'LESIONAL'].copy()
    plot_disease(adata=spatial_adata, obs_name='DISEASE', title='Lesion_skin', save_folder=save_folder)
    # Save coordinates and annotation to .xlsx
    df = pd.DataFrame.from_dict({'UMAP1': spatial_adata.obsm['X_umap'][:, 0],
                                 'UMAP2': spatial_adata.obsm['X_umap'][:, 1],
                                 'DISEASE': spatial_adata.obs['DISEASE'].values})
    df.to_excel(os.path.join(save_folder, 'Plot_infos.xlsx'))


if __name__ == '__main__':
    today = date.today()
    # create saving folder
    output_path = os.path.join("..", "..", "..", "output", "Figure_3A", str(today))
    os.makedirs(output_path, exist_ok=True)

    # Load data:
    pp_st_adata = sc.read(
        os.path.join("..", "..", "..", 'adata_storage', '2022-04-08', 'st_QC_normed_BC_project_PsoADLP.h5'))

    main(save_folder=output_path, spatial_adata=pp_st_adata)
