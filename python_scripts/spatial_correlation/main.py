#!/usr/bin/env python
"""Main script for calculating the (Weighted) Spatial Pearson Correlation for different methods
    File name: main.py
    Author: Christina Hillig
    Credits: Christina Hillig, Ali Farnoud-Niedermayr, Michael Menden
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

from python_scripts.spatial_correlation import density_clustering
from python_scripts.utils import gene_lists

import os
from datetime import date
import scanpy as sc
import anndata


def main(save_folder: str, adata: anndata.AnnData, radius: [int, list], cond_genes: list, genes_resps: dict,
         find_associated_genes: bool = False, corr_method: str = 'pearson', get_plots: bool = False):
    """ Call conditional-based density clustering function

    Parameters
    ----------
    save_folder : str
        path to output directry
    adata : anndata.AnnData
        adata object containing either raw counts or both raw and normed counts.
        Later is necessary if find_associated_genes is True
    radius : int, list
        radius/distance from center spot to surrounding nearest neighbor spots
        e.g.  1 or list [1, 2, 3, ..]
    cond_genes : list
        name of genes to investigate in spots
    genes_resps : dict
        name of cond_genes associated genes such as responder
    find_associated_genes : bool
        if you want to find other than known associated genes of cond_genes
        if True, please provided filtered, normed adata object with sizefactors in .obs as input
    corr_method : str
        which weighted correlation method to use: spearman or pearson (default)
    get_plots : bool
        create evaluation plots

    Returns
    -------

    """
    # parameter
    tissue_layers = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']
    epidermis_layers = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']

    # Run conditional density clustering
    adata, counts_dict, df_stats_responders_in_vs_outlesion_sdc, \
    df_stats_cytokines_responders_in_sdc, df_radius_vs_spearman = density_clustering.main(
        adata=adata, save_folder=save_folder, tissue_types=tissue_layers, epidermis_layers=epidermis_layers,
        radii=radius, get_plots=get_plots, corr_method=corr_method, conditional_genes=cond_genes,
        conditionalgenes_responders=genes_resps, find_associated_genes=find_associated_genes)

    return adata, counts_dict, df_stats_responders_in_vs_outlesion_sdc, \
           df_stats_cytokines_responders_in_sdc, df_radius_vs_spearman


if __name__ == '__main__':
    today = date.today()
    # replace os.environ['PYTHONPATH'].split(os.pathsep)[0] with sys.path[2] -> can be run then in terminal
    path = os.path.join("..", "..")
    # create saving folder in current project path
    savepath = os.path.join(path, "output", "Fig5F-H_Spatial_weighted_correlation", str(today))
    os.makedirs(savepath, exist_ok=True)

    # Load Raw anndata --> used for publication figure 4
    unpp_st_adata = sc.read(
        os.path.join(path, "adata_storage", "2022-04-08", "Spatial Transcriptomics_unpp_cleaned_PsoADLP.h5"))

    # 1. Get cytokines and responders
    conditional_genes, _, conditionalgenes_responders = gene_lists.get_publication_cyto_resps()

    main(save_folder=savepath, adata=unpp_st_adata, cond_genes=conditional_genes,
         genes_resps=conditionalgenes_responders, radius=1, get_plots=True, find_associated_genes=False)
