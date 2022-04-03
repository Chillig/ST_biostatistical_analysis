#! /usr/bin/python
"""Plot UMAP of whole Spatial Transcriptomics data set together with IL-17A positive spots
    File name: UMAP_ST_tissuelayers_IL17A.py
    Author: Christina Hillig
    Date created: December/xx/2020
    Date last modified: May/03/2021
    Python Version: 3.7
"""
from python_scripts.utils import gene_lists, add_observables, get_condition_spots

#import libraries
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


def plot_tissuerlayers_cyto(adata, obs_name, title, save_folder, regions, gene_colors=None):
    """Plot UMAP embedded tissue layers with IL17A on top

    Parameters
    ----------
    adata : annData
    obs_name : str
    title : str
    save_folder : str
    regions : str
    gene_colors: None, list

    Returns
    -------

    """
    if not gene_colors:
        gene_colors = []
        unique_genes = np.unique(adata.obs[obs_name])
        list(unique_genes).remove('Others')
        for gene in signatures.keys():
            if gene in unique_genes:
                gene_colors.append(signatures[gene])

    cyto_adata = adata[adata.obs[obs_name] != "Others"]

    fig, ax = plt.subplots(figsize=figure_size)
    sc.pl.umap(adata, color=regions, use_raw=True, ax=ax, wspace=0.4, show=False,
               size=20, frameon=True, facecolor='white', title="", palette='Set3')
    sc.pl.umap(cyto_adata, color=obs_name, use_raw=True, ax=ax, wspace=0.4, show=False,
               size=50, frameon=True, title="", facecolor='white', vmax=120,
               palette=gene_colors)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    plt.savefig(os.path.join(save_folder, "_".join(['UMAP', title, "Tissuelayers", file_format])))
    plt.close()


def include_cytokine_dp(adata, cytokines, label, save_folder, key, paper_figure):
    """Include double cytokine positive cells in DGE analysis

    Parameters
    ----------
    adata : annData
    cytokines : list
    label : str
    save_folder : str
    key : str
    paper_figure : str

    Returns
    -------

    """
    for cyto in cytokines:
        if "_".join(['cytokine', cyto]) in adata.obs_keys():
            get_condition_spots.get_spots_per_condition_multiple(
                adata=adata, observable="_".join(["cytokine", cyto]), cell_label=label, save_folder=save_folder,
                paper_figure=paper_figure)
            # get_condition_spots.get_spots_per_condition(
            #     adata=adata, observable="_".join(["cytokine", cyto]), save_folder=save_folder, key=key,
            #     paper_figure=paper_figure, cell_label=label)


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


def main(save_folder, spatial_adata):
    """Read out data for ST DGE Analysis and create UMAPs for Figure 3A

    :return:
    """
    spatial_cluster_label = 'tissue_type'

    # load data
    cytokines, allinone, cytoresps_dict = gene_lists.get_publication_cyto_resps()

    # remove all spots without a tissue label
    spatial_adata = spatial_adata[spatial_adata.obs[spatial_cluster_label] != 'Unknown']

    # 1. get observable for cytokine genes
    spatial_adata, obs_name = add_observables.convert_variable_to_observable(
        adata=spatial_adata, gene_names=cytokines, task='cell_gene', label='celltype', condition=None)

    # 2. Highlight tissues epidermis and dermis + cytokines and for each single cytokine
    plot_tissuerlayers_cyto(adata=spatial_adata, obs_name='cytokine_IL17A', title='Wholedataset_IL17A',
                            save_folder=save_folder, regions=spatial_cluster_label)


if __name__ == '__main__':
    today = date.today()
    # create saving folder
    output_path = os.path.join("..", "..", "..", "output", "Figure_1B", str(today))
    os.makedirs(output_path, exist_ok=True)

    # Load data:
    pp_st_adata = sc.read(os.path.join("..", "..", "..", 'adata_storage', '2020-12-04_Visium_Data_QC_BC_clustered.h5'))

    main(save_folder=output_path, spatial_adata=pp_st_adata)
