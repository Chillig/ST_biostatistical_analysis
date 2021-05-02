#!/usr/bin/env python
"""Plot single cell data clusters
    File name: SuppFig4A__SC_UMAP_clusters.py
    Author: Christina Hillig
    Date created: March/xx/2021
    Date last modified: April/30/2021
    Python Version: 3.7
"""
from python_scripts.utils import gene_lists, add_observables

import scanpy as sc
import os
from datetime import date
from collections import OrderedDict
import numpy as np
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

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
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "_".join(['UMAP_annotated_clusters', key, title, ".pdf"])))
    plt.close()


def apply_clusteralgo(adata, algorithm, resolution):
    key = "".join([algorithm, '_r', str(resolution)])
    # clustering
    if "louvain" in algorithm:
        sc.tl.louvain(adata, key_added=key, resolution=resolution, random_state=10)
    else:
        sc.tl.leiden(adata, key_added=key, resolution=resolution)

    # evaluate using silhouette scores
    # Score ranges from −1 to +1
    # -> high value indicates that object is well matched to its own cluster and poorly matched to neighboring clusters
    print(silhouette_score(adata.obsm['X_pca'], adata.obs[key]))
    # silhouette_score(adata.uns['neighbors']['distances'].todense(), adata.obs[key], metric='precomputed')

    return adata, key


def annotate_cluster(adata, cluster_algorithm, resolution):
    """
    0: lymphocytes
    1: antigen presenting cells
    2: fibroblasts
    3: smooth muscle cells
    4: granulocytes
    5: keratinocytes

    :param adata:
    :param cluster_algorithm:
    :param resolution:
    :return:
    """
    obs_key = "".join([cluster_algorithm, '_r', str(resolution)])

    manual_clusterlabels = OrderedDict()
    manual_clusterlabels['0'] = 'Lymphocytes'
    manual_clusterlabels['1'] = 'Antigen presenting cells'
    manual_clusterlabels['2'] = 'Fibroblasts'
    manual_clusterlabels['3'] = 'Smooth muscle cells'
    manual_clusterlabels['4'] = 'Granulocytes'
    manual_clusterlabels['5'] = 'Keratinocytes'

    adata.obs["cluster_labels"] = 'label'
    clusters = np.unique(adata.obs[obs_key])
    for cluster in clusters:
        m_cluster = adata.obs[obs_key] == cluster
        adata.obs['cluster_labels'][m_cluster] = manual_clusterlabels[cluster]

    return adata


def add_tissue_obs(adata):
    # 1 Select tissue types of interest
    tissue_types = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS',
                    'DERdepth1', 'DERdepth2', 'DERdepth3', 'DERdepth4', 'DERdepth5', 'DERdepth6', 'DERdepth7',
                    'INTERFACE']

    adata.obs['tissue_type'] = 'Unknown'
    adata.obs['tissue_type'] = adata.obs['tissue_type'].astype('<U16')
    for tissue in tissue_types:
        m_tissue = adata.obs[tissue] == 1
        adata.obs['tissue_type'][m_tissue] = tissue

    return adata


def add_disease_healthy_obs(adata):
    # 1 Select tissue types of interest
    diagnosis = ['PSO', 'AE', 'LICHEN', 'PRP']

    adata.obs['healthy_disease'] = 'NON LESIONAL'
    adata.obs['healthy_disease'] = adata.obs['healthy_disease'].astype('<U16')
    for disease in diagnosis:
        m_disease = (adata.obs[disease] == 1) & (adata.obs['NON LESIONAL'] == 0)
        adata.obs['healthy_disease'][m_disease] = disease

    return adata


def main(save_folder, pp_adata, cluster_algorithm):
    """
    1. scRNAseq data set
    Read ou pre-processed Count matrix and apply Leiden clustering with resolution r = 0.1 for scRNAseq data set
    Annotate clusters manually

    2. ST data set
    Visualise count matrix with tissue types

    :return:
    """

    # 1. load gene list
    cytokines, allinone, cytoresps_dict = gene_lists.get_publication_cyto_resps()

    # 2. Get observable for cytokine genes
    pp_adata, _ = add_observables.convert_variable_to_observable(
        adata=pp_adata, gene_names=cytokines, task='cell_gene', label='celltype', condition=None)

    # 3. Apply cluster algorithm
    pp_adata, key = apply_clusteralgo(adata=pp_adata, algorithm=cluster_algorithm, resolution=0.1)

    # 4. Annotate clusters with expert opinion - before the best resolution r=0.1 was identified
    pp_adata = annotate_cluster(adata=pp_adata, cluster_algorithm=cluster_algorithm, resolution=0.1)

    # 5. Plot UMAP scRNAseq data
    visualise_clusters(adata=pp_adata, save_folder=save_folder, key='cluster_labels', title="SC")


if __name__ == '__main__':
    today = date.today()
    # create saving folder
    output_path = os.path.join("..", "..", "..", "output", "SupplFigure_4A", str(today))
    os.makedirs(output_path, exist_ok=True)

    # Load data:
    # Use merged scRNAseq samples for publication
    pp_adata_sc = sc.read(os.path.join("..", "..", "..", 'adata_storage', '2020-10-19',
                                       'sc_adata_minumi_600_maxumi_25000_mg_500_mc_20_mt_0_25.h5'))

    unsupvised_cluster_algo = 'leiden'

    main(save_folder=output_path, pp_adata=pp_adata_sc, cluster_algorithm=unsupvised_cluster_algo)
