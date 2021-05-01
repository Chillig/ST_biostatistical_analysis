from python_scripts.utils import gene_lists, add_observables, get_condition_spots

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
            get_condition_spots.get_spots_per_condition(
                adata=adata, observable="_".join(["cytokine", cyto]), save_folder=save_folder, key=key,
                paper_figure=paper_figure, cell_label=label)


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
    """

    :param adata:
    :param tissue_label:
    :return:
    """
    m_epidermis = np.array(
        adata.obs[tissue_label].values)[:, np.newaxis] == np.array(
        ['upper EPIDERMIS', 'basal EPIDERMIS', 'middle EPIDERMIS', 'INTERFACE'])[np.newaxis, :]
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
    """Read out data for ST and scRNA-seq DGE Analysis and create UMAPs for Figure 3A/E and Suppl. Figures 3

    :return:
    """
    spatial_cluster_label = 'tissue_type'

    # load data
    cytokines, allinone, cytoresps_dict = gene_lists.get_publication_cyto_resps()
    leukocyte_markers = gene_lists.leukocyte_markers()

    # remove all spots without a tissue label
    spatial_adata = spatial_adata[spatial_adata.obs[spatial_cluster_label] != 'Unknown']

    # 1. get observable for cytokine genes
    spatial_adata, obs_name = add_observables.convert_variable_to_observable(
        adata=spatial_adata, gene_names=cytokines, task='cell_gene', label='celltype', condition=None)

    spatial_adata, _ = add_observables.convert_variable_to_observable(
        adata=spatial_adata, gene_names=leukocyte_markers,
        task='cell_gene', label='celltype', condition=None)

    # # 2. Read out counts and metaData for DGE Analysis including double positive cytokine cells
    # 2.1 Read out only leukocytes spots by 'CD2', 'CD3D', 'CD3E', 'CD3G', 'CD247' and 'PTPRC' surface markers
    adata_leukocytes = get_celltypes_data(spatial_adata, genes=leukocyte_markers)

    # 2.2 Merge layers of epidermis and save it as epidermis and merge dermis depths and save it as dermis
    adata_leukocytes = get_tissueregions(adata=adata_leukocytes, tissue_label=spatial_cluster_label)

    # 3. Highlicht tissues epidermis and dermis + cytokines and for each single cytokine
    plot_tissueregions_cyto(adata=adata_leukocytes, obs_name='cytokine_IL13', title='Leukocytes_IL13',
                            save_folder=save_folder)
    plot_tissueregions_cyto(adata=adata_leukocytes, obs_name='cytokine_IFNG', title='Leukocytes_IFNG',
                            save_folder=save_folder)

    # 4. Read out all leukocyte psotive spots
    include_cytokine_dp(adata=adata_leukocytes, cytokines=cytokines, save_folder=save_folder,
                        label=spatial_cluster_label, key='ST', paper_figure='3AC_Leukocytes')


if __name__ == '__main__':
    today = date.today()
    wd_path = os.environ['PYTHONPATH'].split(os.pathsep)[0]
    # create saving folder
    output_path = os.path.join(wd_path, "output", "SupplFigure_3AC", str(today))
    os.makedirs(output_path, exist_ok=True)

    # Load data:
    pp_st_adata = sc.read(os.path.join(wd_path, 'adata_storage/2020-12-04_Visium_Data_QC_BC_clustered.h5'))

    main(save_folder=output_path, spatial_adata=pp_st_adata)
