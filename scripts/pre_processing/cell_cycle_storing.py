import scripts.pre_processing.plot_functions.plots_preprocessing as prep_plt

import scanpy as sc
import pandas as pd
import numpy as np


def score_cell_cycle(cc_genes_file, adata, save_folder, raw=False):
    """Identify in which cell cycle a single cells is using cell cycle genes from  Macosko et al., Cell 161 (2015)
    This function should be only applied on single cell resolved data and not on (mini-)bulk RNAseq data

    Parameters
    ----------
    cc_genes_file : str
        path to cc_genes_file
    adata : annData
    save_folder : str
    raw : bool
        if raw annData is provided

    Returns
    -------

    """
    # load gene list from Macosko et al., Cell 161 (2015)
    cc_genes = pd.read_table(cc_genes_file, delimiter=';')

    # drops out S and G2.M cell cycle phases
    if 'human_Cell_cycle_genes.csv' in cc_genes_file:
        # For human_cell_cyle.csv file:
        s_genes = cc_genes['geneID'].values[cc_genes['phase'].values == 'S']
        g2m_genes = cc_genes['geneID'].values[cc_genes['phase'].values == 'G2/M']

        # Get gene names by finding gene_id (ENSG...) in adata variables
        s_genes_hs_ens = adata.var_names[np.isin(adata.var['gene_id'].values, s_genes)]
        g2m_genes_hs_ens = adata.var_names[np.isin(adata.var['gene_id'].values, g2m_genes)]
    else:
        s_genes = cc_genes['S'].dropna()
        g2m_genes = cc_genes['G2.M'].dropna()

        # Get gene names by finding gene name in adata variables
        s_genes_hs_ens = adata.var_names[np.isin(adata.var_names, s_genes)]
        g2m_genes_hs_ens = adata.var_names[np.isin(adata.var_names, g2m_genes)]

    # save genes in adata
    sc.tl.score_genes_cell_cycle(adata=adata, s_genes=s_genes_hs_ens, g2m_genes=g2m_genes_hs_ens)

    # Plot phases
    prep_plt.plot_cell_cycle_cluster(adata=adata, save_folder=save_folder, raw=raw)

    return adata
