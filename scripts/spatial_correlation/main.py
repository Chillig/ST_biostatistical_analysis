#!/usr/bin/env python
"""Main script for calculating the (Weighted) Spatial Pearson Correlation for different methods
    File name: main.py
    Author: Christina Hillig
    Credits: Christina Hillig, Ali Farnoud-Niedermayr, Michael Menden
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

from scripts.spatial_correlation import density_clustering

import os
from datetime import date
import scanpy as sc

if __name__ == '__main__':
    today = date.today()

    # parameter
    tissue_layers = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']
    radius = 1  # or list [1, 2, 3, ..]
    get_plots = False

    # create saving folder in current project path
    savepath = os.path.join(os.environ['PYTHONPATH'].split(os.pathsep)[0], "results", "Correlation", str(today))
    os.makedirs(savepath, exist_ok=True)

    # Load Raw anndata --> used for publication figure 4
    adata_path = os.environ['PYTHONPATH'].split(os.pathsep)[0]
    adata = sc.read(os.path.join(adata_path, "adata_storage/2020-10-06/adata_P15509_P16357_wo_4_7_unpp.h5"))

    density_clustering.main(adata=adata, save_folder=savepath, tissue_types=tissue_layers,
                            radii=radius, get_plots=get_plots)
