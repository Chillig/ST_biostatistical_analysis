import python_scripts.pre_processing.plot_functions.plots_preprocessing as pp_plt

import scanpy as sc
import numpy as np


def find_highly_variable_genes(adata, type_dataset, num_top_genes, save_folder, raw=False):
    """Find highly variable genes in data set which are the most informative ones
    --> information of underlying biological variation
    Parameters
    ----------
    adata : annData
    type_dataset : str
    num_top_genes : int
    save_folder : str
    raw : bool

    Returns
    -------
    adata : annData

    """
    sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=num_top_genes, batch_key='specimen')
    print('\n', 'Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))

    pp_plt.plot_highly_variable_genes(adata, type_dataset, save_folder, raw=raw)

    # plot UMI-count of Housekeeping gene such as GAPDH to check if distribution follows an equal distribution
    pp_plt.plot_hkg_dist(adata=adata, save_folder=save_folder)
    pp_plt.plot_hvg_hkg(adata=adata, save_folder=save_folder, highly_variable_genes=True)

    return adata
