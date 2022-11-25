#!/usr/bin/env python
"""Plot single cell data clusters
    File name: SuppFig4B__SC_UMAP_IFNG.py
    Author: Christina Hillig
    Date created: March/xx/2021
    Date last modified: April/30/2021
    Python Version: 3.7
"""
from python_scripts.utils import gene_lists, add_observables, get_condition_spots

from datetime import date
import scanpy as sc
import numpy as np
import pandas as pd
import os
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


def add_text(adata, obsname, ax, xpos=0.6, max_ypos=0.95):
    cells = np.unique(adata.obs[obsname])
    for ind, cell in enumerate(cells):
        numcells = adata[adata.obs[obsname] == cell].shape[0]
        ax.annotate("No. {} cells: {}".format(cell, numcells),
                    xy=(xpos, max_ypos - (ind * 0.03)), xycoords='axes fraction', fontsize=text_fontsize)


def plot_annotated_cells(adata, color, paper_figure, save_folder,
                         key='SC', title=None, xpos=0.6, ypos=0.95, palette=None):
    if "Spatial" in key:
        size_others = 20
        size_color = 50
    else:
        size_others = 120000 / adata.shape[0]
        size_color = 120000 / adata.shape[0]

    # Remove Others from adata object
    sup_adata = adata[adata.obs[color] == "Others"]
    cyto_adata = adata[adata.obs[color] != "Others"]

    if not palette:
        if len(np.unique(adata.obs[color])) > 2:
            palette = list(signatures.values())
        else:
            gene_name = color.split("_")
            if gene_name != 'others':
                if gene_name[0] in signatures:
                    palette = [signatures[gene_name[0]]]
                else:
                    palette = [signatures[gene_name[-1]]]
            else:
                palette = None

    fig = plt.figure(facecolor='w', edgecolor='k', figsize=figure_size)
    fig.subplots_adjust(bottom=0.025, left=0.025, top=0.975, right=0.975)
    ax = fig.add_subplot(1, 1, 1)
    sc.pl.umap(adata=sup_adata, color=color, use_raw=True, ax=ax, wspace=0.4, size=size_others, frameon=True,
               show=False, title=" ", palette=["lightgrey"])
    sc.pl.umap(cyto_adata, color=color, use_raw=True, ax=ax, wspace=0.4, size=size_color, frameon=True,
               show=False, title=" ", palette=palette)

    # add text of No. cells per cell type
    add_text(adata=adata, obsname=color, ax=ax, xpos=xpos, max_ypos=ypos)
    ax.set_title(paper_figure, loc="left", fontsize=title_fontsize)

    # DOnt show top and right frame lines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    plt.savefig(os.path.join(save_folder, "_".join([key, title, file_format])))
    plt.close()


def include_cytokine_dp(adata, cytokines, label, save_folder, key, paper_figure):
    """
    Include all double cytokine positive cells

    :param adata: [annData]
    :param cytokines: [list]
    :param label: [string]
    :param save_folder: [string]
    :param key:
    :param paper_figure:
    :return:
    """

    for cyto in cytokines:
        if "_".join(['cytokine', cyto]) in adata.obs_keys():
            if 'SC' in key:
                get_condition_spots.get_spots_per_condition(
                    adata=adata, observable="_".join(["cytokine", cyto]), save_folder=save_folder, key=key,
                    paper_figure=paper_figure, cell_label=label)
            else:
                get_condition_spots.get_spots_per_condition_multiple(
                    adata=adata, observable="_".join(["cytokine", cyto]), cell_label=label, save_folder=save_folder,
                    paper_figure=paper_figure)


def exclude_cytokine_dp(adata, cytoresps_dict):
    """Exclusively label double positive cells as double positive

    Parameters
    ----------
    adata : annData
    cytoresps_dict : dict

    Returns
    -------
    adata : annData
    obs_name : str

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
    """Get adata object containing only those cells/spots which express genes of interest

    Parameters
    ----------
    adata : annData
    genes : str, list

    Returns
    -------
    adata : annData

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


def main(save_folder, adata):
    """Read out spots for DGE analysis and create UMAP of single cell RNA-seq data for Suppl. Figures 4B

    Parameters
    ----------
    save_folder : str
    adata : annData

    Returns
    -------

    """
    sc_cluster_obs = 'cluster_labels'

    # 1. load gene list
    cytokines, allinone, cytoresps_dict = gene_lists.get_publication_cyto_resps()
    leukocyte_markers = gene_lists.leukocyte_markers()

    # 2. get observable for cytokine genes
    adata, obs_name = add_observables.convert_variable_to_observable(
        adata=adata, gene_names=cytokines, task='cell_gene', label='celltype', condition=None)

    # Only Leukocytes:
    # 3. Read out counts and metaData for DGE Analysis including double positive cytokine cells
    # 3.1 Read out only T-cell spots by CD2 surface markers
    adata_leukocytes = get_celltypes_data(adata, genes=leukocyte_markers)

    # 3.3 Read out all leukocyte cells
    include_cytokine_dp(adata=adata_leukocytes, cytokines=cytokines, save_folder=save_folder,
                        label=sc_cluster_obs, key='SC_merged', paper_figure='SC')
    # 3.3 Read out all leukocyte cells but exclude double positive cytokine cells
    adata_leukocytes, obs_name = exclude_cytokine_dp(adata=adata_leukocytes, cytoresps_dict=cytoresps_dict)

    # Plot cytokines and highlight double positive
    plot_annotated_cells(adata=adata_leukocytes, color='cytokines_others', paper_figure='',
                         save_folder=save_folder, key='SC', title="Leukocytes_IL17A_IFNG",
                         xpos=0.02, ypos=0.95, palette=["#ff7f00", "#377eb8", 'purple'])

    # Add cytokine label to adata and Plot: Highlight cytokines
    adata_leukocytes = add_observables.add_columns_genes(adata=adata_leukocytes, genes='IFNG', label='IFNG')

    """ Suppl. Figure 6A: Highlight IFN-g """
    plot_annotated_cells(adata=adata_leukocytes, color='IFNG_label', paper_figure='6A', save_folder=save_folder,
                         key='SC', title="Leukocyte_IFNG", xpos=0.02, ypos=0.95)

    # Save coordinates and annotation to .xlsx
    df = pd.DataFrame.from_dict({'UMAP1': adata_leukocytes.obsm['X_umap'][:, 0],
                                 'UMAP2': adata_leukocytes.obsm['X_umap'][:, 1],
                                 'IFNG_label': adata_leukocytes.obs['IFNG_label'].values})
    df.to_excel(os.path.join(save_folder, 'Plot_infos.xlsx'))


if __name__ == '__main__':
    today = date.today()
    # create saving folder
    output_path = os.path.join("..", "..", "..", "output", "SupplFigure_6A", str(today))
    os.makedirs(output_path, exist_ok=True)

    # Load data:
    # Use merged scRNAseq samples for publication
    clustered_adata_sc = sc.read(
        os.path.join("/Users/christina.hillig/PycharmProjects/Cellranger_analysis/Spatial_publication/Input_data",
                     '2020-12-04_SC_Data_QC_clustered.h5'))

    main(save_folder=output_path, adata=clustered_adata_sc)
