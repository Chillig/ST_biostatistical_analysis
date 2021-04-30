#!/usr/bin/env python
"""Read out of cytokine counts in RNA-seq data and visualisation of data in PCA and UMAP plots
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
import seaborn as sns
from sklearn import decomposition
import scanpy as sc

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


def get_counts(data, metadata, cytokine, diagnosis):
    """

    Parameters
    ----------
    data : pandas.Dataframe
    metadata : pandas.Dataframe
    cytokine: str
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

    df = pd.DataFrame(columns=['counts', 'skin'])
    df['counts'] = diagnosis_data.loc[cytokine].values
    df['skin'] = metadata['skin'].values[m_diagnosis]

    return df


def plot_count_distribution(counts, cytokine, save_folder, key, title):
    """Plot count distribution of RNA-seq data

    Parameters
    ----------
    counts : pandas.Dataframe
    cytokine : str
    save_folder : str
    key : str
    title : str

    Returns
    -------

    """
    fig, ax = plt.subplots(figsize=fig_size)
    ax.set_axisbelow(True)
    ax.grid(color='gray', linestyle='dashed')
    # Lesional
    m_lesional = counts['skin'] == 'lesional'
    num_biopsies_lesional = len(counts['counts'].values[m_lesional])
    ax.scatter(x=np.arange(0, 2 * num_biopsies_lesional, 2), y=counts['counts'].values[m_lesional], s=8,
               c='k', label='lesional')
    # Non Lesional
    m_nonlesional = counts['skin'] == 'non-lesional'
    num_biopsies_nonlesional = len(counts['counts'].values[m_nonlesional])
    ax.scatter(x=np.arange(1, 2 * num_biopsies_nonlesional, 2), y=counts['counts'].values[m_nonlesional], s=8,
               c='grey', label='Non-lesional')
    ax.set_xlabel('Biopsy', fontsize=xy_fontsize)
    if cytokine == "IFNG":
        ax.set_ylabel(r'IFN$\gamma$ counts', fontsize=xy_fontsize)
    else:
        ax.set_ylabel(" ".join([cytokine, 'counts']), fontsize=xy_fontsize)
    ax.tick_params(labelsize=xy_ticks)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.annotate("Total No. biopsies: {}".format(int(len(counts.values))), xy=(0.68, 1.05),
                xycoords='axes fraction')
    ax.legend(loc='best', title_fontsize=title_fontsize, fontsize=legend_fontsize)

    plt.savefig(os.path.join(save_folder, "_".join([key, title, fileformat])))
    plt.close()


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

    # logaritmize
    sc.pp.log1p(adata)
    adata.raw = adata
    # TODO apply batch effect correction using pyComBat

    return adata


def plot_scree(pca, save_folder, n_comps):
    pc_values = np.arange(pca.n_components_) + 1

    fig, ax = plt.subplots(figsize=fig_size)
    ax.plot(pc_values, pca.explained_variance_ratio_, 'ro-', linewidth=2)
    ax.set_title('Scree Plot', fontsize=title_fontsize)
    ax.set_xlabel('Principal Component', fontsize=xy_fontsize)
    ax.set_ylabel('Proportion of Variance Explained', fontsize=xy_fontsize)
    plt.savefig(os.path.join(save_folder, "".join([str(n_comps), "Bulk_Elbowplot.png"])))
    plt.close()


def apply_pca(adata, df_counts, save_folder, n_comps=50):
    """

    Parameters
    ----------
    adata : annData
        RNA-seq object
    df_counts : pandas.Dataframe
        input = sub_bulk_adata.X.T
    save_folder : str
        path of output folder
    n_comps : int
        Number of Principal Components

    Returns
    -------

    """
    pca = decomposition.PCA(n_components=n_comps, svd_solver='arpack')
    pca.fit(df_counts.T)
    pca_counts = pca.transform(df_counts.T)

    plot_scree(pca=pca, save_folder=save_folder, n_comps=n_comps)

    fig, ax = plt.subplots(figsize=fig_size)
    ax.plot(np.cumsum(pca.explained_variance_ratio_), zorder=2)
    ax.axhline(y=0.75, c='black', linestyle='dashed', lw=0.5, zorder=1)
    ax.set_xlabel('number of components', fontsize=xy_fontsize)
    ax.set_ylabel('cumulative explained variance', fontsize=xy_fontsize)
    plt.savefig(os.path.join(save_folder, "".join([str(n_comps), "Bulk_Variance.png"])))
    plt.close()

    # Save PC in .obsm of annData object
    adata.obsm['pca'] = pca_counts

    return pca_counts


def apply_scanpy_pca(adata, save_folder, n_comps=50):
    sc.pp.pca(adata, n_comps=n_comps, use_highly_variable=False, svd_solver='arpack')

    sc.pl.pca_variance_ratio(adata, log=False, show=False, save=False)
    plt.savefig(os.path.join(save_folder, "".join([str(n_comps), "Bulk_Variance.png"])))
    plt.close()

    return adata


def plot_pca(adata, observable, save_folder):
    # signatures = OrderedDict()
    # # publication
    # signatures["IFNG"] = "#ff7f00"  # orange LICHEN
    # signatures["IL13"] = "#e41a1c"  # red AE
    # signatures["IL17A"] = "#377eb8"  # blue PSO
    fig, ax = plt.subplots(figsize=fig_size)
    if len(np.unique(adata.obs[observable])) > 4:
        sc.pl.pca(adata, color=observable, ax=ax, show=False)
    elif observable == 'diseases':
        sc.pl.pca(adata, color=observable, ax=ax, show=False,
                   palette=["#e41a1c", 'darkgreen', '#ff7f00', "#377eb8", "mediumblue"])
    else:
        sc.pl.pca(adata, color=observable, ax=ax, show=False, palette=["#e41a1c", '#ff7f00', "#377eb8", "mediumblue"])
    ax.set_xlabel("PC1 ({:.0f}% Variance explained)".format(adata.uns["pca"]["variance_ratio"][0] * 100))
    ax.set_ylabel("PC2 ({:.0f}% Variance explained)".format(adata.uns["pca"]["variance_ratio"][1] * 100))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "".join(["PCA-Bulk", observable.capitalize(), fileformat])))
    plt.close()


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
    # Read out IL17A per pso biopsy
    il17a_counts = get_counts(data=bulk_rnaseq, metadata=metadata, cytokine='IL17A', diagnosis='psoriasis')

    # Read out IFNG per LICHEN biopsy
    ifng_counts = get_counts(data=bulk_rnaseq, metadata=metadata, cytokine='IFNG', diagnosis='lichen ruber')

    # Read out IL13 per AE biopsy
    il13_counts = get_counts(data=bulk_rnaseq, metadata=metadata, cytokine='IL13', diagnosis='Atopic eczema')

    # Read out all cytokines per PRP biospy
    il17a_prp_counts = get_counts(data=bulk_rnaseq, metadata=metadata,
                                  cytokine='IL17A', diagnosis='pytiriasis rubra pilaris')
    ifng_prp_counts = get_counts(data=bulk_rnaseq, metadata=metadata,
                                 cytokine='IFNG', diagnosis='pytiriasis rubra pilaris')
    il13_prp_counts = get_counts(data=bulk_rnaseq, metadata=metadata,
                                 cytokine='IL13', diagnosis='pytiriasis rubra pilaris')

    # Save as .csv file
    il17a_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__PSO__IL17A.csv"))
    ifng_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__LICHEN__IFNG.csv"))
    il13_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__AD__IL13.csv"))
    # PRP
    il17a_prp_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__PRP__IL17A.csv"))
    ifng_prp_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__PRP__IFNG.csv"))
    il13_prp_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__PRP__IL13.csv"))

    print("Plot")
    plot_count_distribution(counts=il17a_counts, cytokine='IL17A',
                            save_folder=save_folder, key="Bulk_RNAseq", title="_".join(['IL17A', 'PSO_distribution']))
    plot_count_distribution(counts=ifng_counts, cytokine='IFNG',
                            save_folder=save_folder, key="Bulk_RNAseq", title="_".join(['IFNG', 'LICHEN_distribution']))
    plot_count_distribution(counts=il13_counts, cytokine='IL13',
                            save_folder=save_folder, key="Bulk_RNAseq", title="_".join(['IL13', 'AD_distribution']))
    plot_count_distribution(counts=il17a_prp_counts, cytokine='IL17A', save_folder=save_folder, key="Bulk_RNAseq",
                            title="_".join(['IL17A', 'PRP_distribution']))
    plot_count_distribution(counts=ifng_prp_counts, cytokine='IFNG', save_folder=save_folder, key="Bulk_RNAseq",
                            title="_".join(['IFNG', 'PRP_distribution']))
    plot_count_distribution(counts=il13_prp_counts, cytokine='IL13', save_folder=save_folder, key="Bulk_RNAseq",
                            title="_".join(['IL13', 'PRP_distribution']))

    # plot log2(expr + 1) histogram
    fig_genes = plt.figure(facecolor='w', edgecolor='k', figsize=fig_size)
    fig_genes.subplots_adjust(bottom=0.025, left=0.025, top=0.975, right=0.975)
    ax = fig_genes.add_subplot(1, 1, 1)
    sns.distplot(np.log2(bulk_rnaseq.values + 1), kde=False, bins=50)
    ax.set_ylabel(r'log_{2}(counts + 1)', fontsize=xy_fontsize)
    ax.set_xlabel('Counts', fontsize=xy_fontsize)
    ax.tick_params(labelsize=xy_ticks)
    plt.savefig(os.path.join(save_folder, "Counts_Distribution.pdf"))
    plt.close()


if __name__ == '__main__':
    today = date.today()
    wd_path = os.environ['PYTHONPATH'].split(os.pathsep)[0]
    # create saving folder
    output_path = os.path.join(wd_path, "output", "Figure_2B", str(today))
    os.makedirs(output_path, exist_ok=True)

    # input path
    input_path = "/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/input/bulk_RNAseq"

    # Read bulk-RNAseq count matrix
    bulk_data = pd.read_csv(os.path.join(input_path, "bulkRNA_countMat.txt"), sep='\t')
    # Read in metaData
    meta_data = pd.read_excel(os.path.join(input_path, "bulkRNA_metaData.xlsx"))

    main(save_folder=output_path, bulk_rnaseq=bulk_data, metadata=meta_data)
