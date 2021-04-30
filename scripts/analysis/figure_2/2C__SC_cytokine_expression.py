from scripts.utils import gene_lists, add_observables as ctools

import scanpy as sc
import numpy as np
import pandas as pd
import os
from datetime import date


def get_counts(adata, genes):
    """Get counts of genes

    Parameters
    ----------
    adata : annData
    genes : str, list

    Returns
    -------
    adata: annData

    """
    # get counts of each gene in genes
    for gene in genes:
        obs_counts = "_".join([gene, 'counts'])
        varindex_genes = np.where(adata.var.index[np.newaxis, :] == np.array([gene])[:, np.newaxis])[1]
        if "counts" in adata.layers.keys():
            counts_genes = np.copy(adata.layers['counts'])[:, varindex_genes]
        else:
            counts_genes = np.copy(adata.X)[:, varindex_genes]

        # Add observables to adata
        adata.obs[obs_counts] = counts_genes
    return adata


def get_cd_cytokine_positive_cells(adata, obs_cyto, obs_celltype, celltype):
    """Get counts of CD4/8 surface markers

    Parameters
    ----------
    adata
    obs_cyto
    obs_celltype
    celltype

    Returns
    -------

    """
    obs_celllabel = "_".join(['cytokine', obs_celltype])
    obs_cellcounts = "_".join([obs_celltype, 'counts'])
    # Read out those cells which are not labeled as "Others"
    adata_celltype = adata[adata.obs[obs_celllabel] == celltype].copy()

    # Create mask to identify cytokine positive cells
    m_cyto_positive = adata_celltype.obs[obs_cyto] > 0

    # Sum over those cells from the celltype which are cytokine positive cells
    celltype_cyto_positive = adata_celltype.obs[obs_cellcounts][m_cyto_positive]
    # Sum over cytokine negative cells among the celltype cells
    celltype_cyto_negative = adata_celltype.obs[obs_cellcounts][~m_cyto_positive]

    # Generate additional information for plotting
    total_counts = adata_celltype.obs[obs_cellcounts].sum()
    num_cells = adata_celltype.shape[0]
    data = [celltype_cyto_positive.sum(), celltype_cyto_negative.sum()]
    plot_label = ["".join([obs_cyto.split("_")[0], "+"]), "".join([obs_cyto.split("_")[0], "-"])]

    return data, plot_label, total_counts, num_cells, celltype_cyto_positive, celltype_cyto_negative


def save_cd_counts(adata, cytokine, save_folder):
    """Save counts of CD4 and CD8 surface marker in .csv file

    Parameters
    ----------
    adata : annData
    cytokine : str
    save_folder : str

    Returns
    -------

    """
    obs_cyto = "_".join([cytokine, 'counts'])
    data, plot_label, total_counts, num_cells, cd4_cyto_positive, cd4_cyto_negative = \
        get_cd_cytokine_positive_cells(adata, obs_cyto=obs_cyto, obs_celltype='CD4', celltype='CD4')

    data, plot_label, total_counts, num_cells, cd8a_cyto_positive, cd8a_cyto_negative = \
        get_cd_cytokine_positive_cells(adata, obs_cyto=obs_cyto, obs_celltype='CD8A', celltype='CD8A')

    data, plot_label, total_counts, num_cells, cd8b_cyto_positive, cd8b_cyto_negative = \
        get_cd_cytokine_positive_cells(adata, obs_cyto=obs_cyto, obs_celltype='CD8B', celltype='CD8B')

    # save counts of cell type positive cells to csv file
    cd_adata = adata.copy()[(np.any(adata.obs[['cytokine_CD4', 'cytokine_CD8A', 'cytokine_CD8B']].values != 'Others',
                                    axis=1)) & (adata.obs[obs_cyto] > 0)]
    df = pd.DataFrame(columns=['CD4', 'CD8'], index=cd_adata.obs.index)
    df['CD4'] = cd4_cyto_positive
    df['CD8'] = cd8a_cyto_positive + cd8b_cyto_positive

    df.to_csv(os.path.join(save_folder, "_".join([cytokine, "CD4_CD8.csv"])))


def main(save_folder, adata):
    """
    1. Identify CD4 and CD8 cells
    2. Get counts of CD4 and CD8 cells for those which are also cytokine+ cells
    3. Get counts of CD4 and CD8 cells for those which are cytokine- cells
    4. Show that in a stacked histogram

    Parameters
    ----------
    save_folder : str
    adata : annData

    Returns
    -------

    """
    celltypes = gene_lists.get_exclusive_cd4cd8()

    # 1.1 Identify CD4 and CD8 cells
    for celllabel in celltypes.keys():
        adata, _ = ctools.convert_variable_to_observable(adata=adata, gene_names=celltypes[celllabel], task='cell_gene',
                                                         condition=None,)
    # 1.2 Get counts of cytokine+ cells
    adata = get_counts(adata=adata, genes=['IL17A', 'IFNG', 'IL13'])

    # 2. Illustrate how many cells of the total CD4 or CD8 cells are positive for IL-17:
    # -> to get ratio of CD4/CD8 IL17A positive cells
    # -> Steffi: plot with relative number of IL-17+ cells amongst CD4 and amongst CD8 cells
    save_cd_counts(adata=adata, cytokine='IL17A', save_folder=save_folder)
    save_cd_counts(adata=adata, cytokine='IFNG', save_folder=save_folder)
    save_cd_counts(adata=adata, cytokine='IL13', save_folder=save_folder)


if __name__ == '__main__':
    today = date.today()
    wd_path = os.environ['PYTHONPATH'].split(os.pathsep)[0]
    # create saving folder
    output_path = os.path.join(wd_path, "output", "Figure_2C", str(today))
    os.makedirs(output_path, exist_ok=True)

    # Load data:
    # Use merged scRNAseq samples for publication
    adata_sc = sc.read(os.path.join(wd_path, 'adata_storage/2020-11-30/sc_adata_unpp.h5'))

    main(save_folder=output_path, adata=adata_sc)
