#!/usr/bin/env python
"""Plot UMAP of diagnosis, skin layers, skin layers + cytokine double positive spots
    File name: SuppFig2ABC__ST_UMAP.py
    Author: Christina Hillig
    Date created: December/xx/2020
    Date last modified: May/02/2021
    Python Version: 3.7
"""
from python_scripts.utils import gene_lists, add_observables

import scanpy as sc
import numpy as np
import os
from datetime import date
from collections import OrderedDict
import matplotlib.pyplot as plt


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


def get_color_palette(num_colors):
    if 20 < num_colors < 28:
        palette = sc.pl.palettes.zeileis_28
    elif 28 <= num_colors < 102:
        palette = sc.pl.palettes.godsnot_102
    else:
        palette = sc.pl.palettes.vega_20_scanpy
    return palette


def plot_tissueregions_cyto(adata, obs_name, title, save_folder, gene_colors=None):
    if not gene_colors:
        gene_colors = []
        unique_genes = np.unique(adata.obs[obs_name])
        list(unique_genes).remove('Others')
        for gene in signatures.keys():
            if gene in unique_genes:
                gene_colors.append(signatures[gene])

    cyto_adata = adata[adata.obs[obs_name] != "Others"]

    fig, ax = plt.subplots(figsize=figure_size)
    sc.pl.umap(adata, color='tissue_regions', use_raw=True, ax=ax, wspace=0.4, show=False,
               size=20, frameon=True, facecolor='white', palette=['bisque', 'pink'],
               title="")
    sc.pl.umap(cyto_adata, color=obs_name, use_raw=True, ax=ax, wspace=0.4, show=False,
               size=50, frameon=True, title="", facecolor='white', vmax=120,
               palette=gene_colors)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    plt.savefig(os.path.join(save_folder, "_".join(['UMAP', title, "Skinlayers", file_format])))
    plt.close()


def visualise_clusters(adata, save_folder, key, title):
    """

    :param adata:
    :param save_folder:
    :param key: [string] Name of cluster algorithm and resolution eg: leiden_r1
    :param title:
    :return:
    """
    num_clusters = len(np.unique(adata.obs[key]))
    palette = get_color_palette(num_colors=num_clusters)

    fig = plt.figure(facecolor='w', edgecolor='k', figsize=figure_size)
    ax = fig.add_subplot(1, 1, 1)
    sc.pl.umap(adata, color=key, use_raw=False, palette=palette, frameon=True, show=False, title='', ax=ax)
    # sc.pl.umap(adata, color='Cytokines', use_raw=False, palette=palette, frameon=True, show=False, title='', ax=ax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "_".join(['UMAP_annotated_clusters', key, title, ".pdf"])))
    plt.close()


def exclude_cytokine_dp(adata, cytoresps_dict):
    """Exclusively label double positive cells as double positive

    Parameters
    ----------
    adata : annData
    cytoresps_dict : dict

    Returns
    -------

    """
    cytodict = gene_lists.cyto_asdict()
    # check if all genes are in adata
    cells = list(set(adata.var.index) & set(cytodict.keys()))
    remove_key = np.setdiff1d(list(cytoresps_dict.keys()), cells)
    if len(remove_key) > 0:
        cytodict.pop(remove_key[0], None)
    adata, obs_name = add_observables.convert_variable_to_observable(
        adata=adata, gene_names=cytodict, task='annotate_cells', label='cytokines', condition=[np.all, np.all, np.all])

    return adata, obs_name


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
        varindex_cyto_genes = np.where(adata.var.index == genes)[0]
        counts_cyto = adata.layers["counts"][:, varindex_cyto_genes][:, 0]
    # create mask
    m_cyto = counts_cyto > 0
    m_cyto = np.any(m_cyto, axis=1)
    adata = adata.copy()[m_cyto]

    return adata


def get_tissueregions(adata, tissue_label):
    """

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


def main(save_folder, spatial_adata, spatial_cluster_label: str = 'tissue_layer'):
    """ Read out data for ST and create UMAPs for Suppl. Figures 2E/F/G

    Parameters
    ----------
    save_folder
    spatial_adata
    spatial_cluster_label

    Returns
    -------

    """

    # 1. load gene lists
    cytokines, allinone, cytoresps_dict = gene_lists.get_publication_cyto_resps()
    leukocyte_markers = gene_lists.leukocyte_markers()

    # 2. Read out only lesional spots
    spatial_adata = spatial_adata[spatial_adata.obs['biopsy_type'] == 'LESIONAL'].copy()

    # 3. get observable for cytokine genes and leukocyte markers
    spatial_adata, obs_name = add_observables.convert_variable_to_observable(
        adata=spatial_adata, gene_names=cytokines, task='cell_gene', label='celltype', condition=None)

    spatial_adata, obs_name = add_observables.convert_variable_to_observable(
        adata=spatial_adata, gene_names=leukocyte_markers, task='cell_gene', label='celltype', condition=None)

    # 4. Read out only leukocytes spots by 'CD2', 'CD3D', 'CD3E', 'CD3G', 'CD247' and 'PTPRC' surface markers
    adata_leukocytes = get_celltypes_data(spatial_adata, genes=leukocyte_markers)

    # 5. add observable healthy_disease
    spatial_adata = add_observables.add_disease_healthy_obs(spatial_adata)

    # keys: 'patient', 'biopsy_type', 'disease', 'tissue_type'
    # Suppl Figure 2A
    visualise_clusters(adata=spatial_adata, save_folder=save_folder, key='healthy_disease', title="LESIONAL_Diagnoses")
    # Suppl. Figure 2C
    visualise_clusters(adata=adata_leukocytes, save_folder=save_folder, key='tissue_layer',
                       title="LESIONAL_Leukocytes_tissuelayers")
    # Suppl. Figure 2D
    visualise_clusters(adata=adata_leukocytes, save_folder=save_folder, key='spot_type',
                       title="LESIONAL_Leukocytes_tissuelayers")

    # 6. Read out spots which either have IL17A, IL13 or INFG genes
    adata_leukocytes, obs_name = exclude_cytokine_dp(adata=adata_leukocytes, cytoresps_dict=cytoresps_dict)

    # 7. Merge layers of epidermis and save it as epidermis and merge dermis depths and save it as dermis
    adata_leukocytes = get_tissueregions(adata=adata_leukocytes, tissue_label=spatial_cluster_label)

    # Suppl. Figure 2B
    plot_tissueregions_cyto(adata=adata_leukocytes, obs_name=obs_name, title='LESIONAL_Leukocytes_Cytokines',
                            save_folder=save_folder, gene_colors=["#ff7f00", "#e41a1c", 'darkgoldenrod', 'purple',
                                                                  "#377eb8", 'deeppink'])


if __name__ == '__main__':
    today = date.today()
    # create saving folder
    output_path = os.path.join("..", "..", "..", "output", "SupplFigure_2EFG", str(today))
    os.makedirs(output_path, exist_ok=True)

    # Load data:
    pp_st_adata = sc.read(os.path.join("..", "..", "..", 'adata_storage', '2020-12-04_Visium_Data_QC_BC_clustered.h5'))

    main(save_folder=output_path, spatial_adata=pp_st_adata)
