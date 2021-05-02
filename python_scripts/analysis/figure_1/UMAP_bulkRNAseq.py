#!/usr/bin/env python
"""Read out of cytokine counts in RNA-seq data and visualisation of data in UMAP plots
    File name: bulk_RNAseq_data.py
    Author: Christina Hillig
    Date created: 3/11/2020
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

from datetime import date
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import anndata
import scanpy as sc
import sys

fig_size = (8, 8)
xy_fontsize = 16
xy_ticks = 12
title_fontsize = 18
legend_fontsize = 14
text_fontsize = 14
fileformat = '.pdf'


def subset_data(data, metadata, diagnosis):
    """Subset data to specific diagnosis

    Parameters
    ----------
    data : pandas.Dataframe
    metadata : pandas.Dataframe
    diagnosis : list or str

    Returns
    -------

    """

    if isinstance(diagnosis, list):
        m_diagnosis = np.where(np.array(metadata['diagnosis'])[:, np.newaxis] == np.array(diagnosis)[np.newaxis, :])[0]
    else:
        m_diagnosis = metadata['diagnosis'] == diagnosis
    biopsies_names = data.columns[m_diagnosis]
    diagnosis_data = data[biopsies_names]

    return diagnosis_data, metadata.iloc[m_diagnosis, :]


def filter_bulkdata(adata):
    adata.obs['n_counts'] = adata.X.sum(1)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    adata.obs['n_genes'] = (adata.X > 0).sum(1)

    mt_gene_mask = [gene.startswith('MT-') for gene in adata.var_names]
    adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1) / adata.obs['n_counts']

    print('Total number of cells: {:d}'.format(adata.n_obs))

    sc.pp.filter_cells(adata, min_counts=500)
    print('Number of cells after min count filter: {:d}'.format(adata.n_obs))
    # Threshold for MT-fraction is 20-25%
    adata = adata[adata.obs['mt_frac'] < 0.25]
    print('Number of cells after MT filter: {:d}'.format(adata.n_obs))
    # Min 20 cells - filters out 0 count genes
    sc.pp.filter_genes(adata, min_cells=20)
    print('Number of genes after cell filter: {:d}'.format(adata.n_vars))

    return adata


def appply_preprocessing(adata):
    # Filter and remove low expressed genes
    adata = filter_bulkdata(adata=adata)
    # normalize to depth 10 000
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

    # log-transform
    sc.pp.log1p(adata)
    adata.raw = adata

    return adata


def apply_scanpy_pca(adata, save_folder, n_comps=50):
    sc.pp.pca(adata, n_comps=n_comps, use_highly_variable=False, svd_solver='arpack')

    sc.pl.pca_variance_ratio(adata, log=False, show=False, save=False)
    plt.savefig(os.path.join(save_folder, "".join([str(n_comps), "Bulk_Variance.png"])))
    plt.close()

    return adata


def plot_umap(adata, observable, save_folder):
    fig, ax = plt.subplots(figsize=fig_size)
    if len(np.unique(adata.obs[observable])) > 5:
        sc.pl.umap(adata, color=observable, ax=ax, show=False)
    elif observable == 'diseases':
        sc.pl.umap(adata, color=observable, ax=ax, show=False,
                   palette=["#e41a1c", 'darkgreen', '#ff7f00', "#377eb8", "mediumblue"])
    else:
        sc.pl.umap(adata, color=observable, ax=ax, show=False, palette=["#e41a1c", '#ff7f00', "#377eb8", "mediumblue"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "".join(["UMAP-Bulk", observable.capitalize(), fileformat])))
    plt.close()


def main(save_folder, bulk_rnaseq, metadata):
    """Read out counts of cytokines in RNA-seq data and create PCA and UMAP plots

    Parameters
    ----------
    save_folder : str
    bulk_rnaseq : pandas.Dataframe
    metadata : pandas.Dataframe

    Returns
    -------

    """

    # put sub-bulk RNA-seq into an annData object;  all diagnosis list(np.unique(metadata['diagnosis']))
    sub_bulkrnaseq, sub_metadata = subset_data(
        data=bulk_rnaseq, metadata=metadata,
        diagnosis=['psoriasis', 'lichen ruber', 'Atopic eczema'])
    data_var = pd.DataFrame()
    data_var['gene_name'] = list(sub_bulkrnaseq.index)
    sub_bulk_adata = anndata.AnnData(sub_bulkrnaseq.values.T, obs=sub_metadata, var=data_var, dtype='int64')

    # Combine diagnosis with skin (lesional and non-lesional information) covariate
    m_lesional = sub_bulk_adata.obs['skin'].values == 'lesional'
    sub_bulk_adata.obs['diseases'] = 'healthy'
    sub_bulk_adata.obs['diseases'][m_lesional] = sub_bulk_adata.obs['diagnosis'][m_lesional]
    sub_bulk_adata.obs['diseases'] = sub_bulk_adata.obs['diseases'].astype('category')

    # Apply Pre-processing
    sub_bulk_adata = appply_preprocessing(adata=sub_bulk_adata)

    # Apply PCA
    sub_bulk_adata = apply_scanpy_pca(adata=sub_bulk_adata, save_folder=save_folder, n_comps=50)
    n_comps = 6  # Determined via explained variance plot using the elbow method
    sub_bulk_adata = apply_scanpy_pca(adata=sub_bulk_adata, save_folder=save_folder, n_comps=n_comps)

    try:
        sc.pp.neighbors(sub_bulk_adata, n_pcs=n_comps, n_neighbors=15, knn=True)
    except AssertionError:
        sc.pp.neighbors(sub_bulk_adata, n_pcs=n_comps, n_neighbors=15, knn=True)
    sc.tl.umap(sub_bulk_adata)

    plot_umap(adata=sub_bulk_adata, observable='diagnosis', save_folder=save_folder)
    plot_umap(adata=sub_bulk_adata, observable='diseases', save_folder=save_folder)
    plot_umap(adata=sub_bulk_adata, observable='age', save_folder=save_folder)
    plot_umap(adata=sub_bulk_adata, observable='gender', save_folder=save_folder)
    plot_umap(adata=sub_bulk_adata, observable='skin', save_folder=save_folder)


if __name__ == '__main__':
    today = date.today()
    # create saving folder
    output_path = os.path.join("..", "..", "..", "output", "Figure_1_Bulk_RNAseq", str(today))
    os.makedirs(output_path, exist_ok=True)

    # input path
    input_path = os.path.join("..", "..", "..", "input", "bulk_RNAseq")

    # Read bulk-RNAseq count matrix
    bulk_data = pd.read_csv(os.path.join(input_path, "bulkRNA_countMat.txt"), sep='\t')
    # Read in metaData
    meta_data = pd.read_excel(os.path.join(input_path, "bulkRNA_metaData.xlsx"))

    main(save_folder=output_path, bulk_rnaseq=bulk_data, metadata=meta_data)
