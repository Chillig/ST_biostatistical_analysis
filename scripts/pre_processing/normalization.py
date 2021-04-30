from scripts.pre_processing.plot_functions import plots_preprocessing
from scripts.utils import helper_tools as ht

import scanpy as sc
import numpy as np

import rpy2.robjects.packages as rpackages
# import rpy2.rinterface_lib as rinfacelib

from rpy2.robjects import numpy2ri
from rpy2.robjects import pandas2ri

numpy2ri.activate()
pandas2ri.activate()


def first_est_sizefactors(adata, save_folder, counts_per_cell_after=1e6, n_comps=15, resolution=0.5, raw=False):
    """Estimate size factors for normalisation of the data

    Normalization method implemented in scran package performs best.
    Method requires a coarse clustering input to improve size factor estimation performance.
    ==> simple pre-processing approach and cluster the data at a low resolution to get input for size factor estimation

    Parameters
    ----------
    adata : annData
    save_folder : str
    counts_per_cell_after : int
        Number of UMI-counts after normalisation in each cell or spot
        value is used as a general scaling of the data set
    n_comps : int
        Number of principle components (PC) to keep (default: 15)
    resolution : float
        low resolution for cluster algorithm to get an input for size factor estimation
    raw : bool
        if raw annData is provided

    Returns
    -------
    adata : annData

    """
    # import R libraries:
    scranlib = rpackages.importr('scran')

    adata_pp = adata.copy()

    # Use annotations from pathologist instead of clusters
    obs_keys = list(adata.obs_keys())
    # get all manual annotations by extracting all keys with upper case characters
    annotations = [char for char in obs_keys if any(c.isupper() for c in char)]

    # check if manual annotations in data set
    if len(annotations) > 0:
        # save manual annotations in one array as groups
        adata_pp = ht.store_categories_as_clusters(adata_pp, annotations)
    else:
        sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=counts_per_cell_after)
        sc.pp.log1p(adata_pp)
        sc.pp.pca(adata_pp, n_comps=n_comps)
        sc.pp.neighbors(adata_pp)
        sc.tl.louvain(adata_pp, key_added='groups', resolution=resolution)

    # Preprocess variables for scran normalization
    input_groups = adata_pp.obs['groups']
    data_mat = adata.X.T

    # estimateSizeFactors: This calculates the relative library depth of each sample or group (cell types)?
    size_factors = scranlib.computeSumFactors(data_mat, clusters=input_groups)

    del adata_pp

    # Visualize the estimated size factors
    adata.obs['size_factors'] = size_factors

    plots_preprocessing.plot_sc_scatter(adata=adata, save_folder=save_folder,
                                        obs=["size_factors", "n_counts", "n_genes"],
                                        ax_xlabel="size factors", ax_ylabel=["No. counts", "No. genes"],
                                        n_rows=1, n_cols=2, raw=raw)

    plots_preprocessing.plot_sc_dist(adata=adata, save_folder=save_folder, observable="size_factors",
                                     bins=50, kde=False,
                                     title="Distribution of size factors",
                                     ax_xlabel="size factor", ax_ylabel="counts", raw=raw)

    return adata


def normalize_by_sizefactor(adata, adata_path_filenames):
    """Normalising the count matrix using sizefactors
    The basic preprocessing includes assuming all size factors are equal
    (library size normalization to counts per million - CPM) and log-transforming the count data

    Parameters
    ----------
    adata : annData
    adata_path_filenames : str

    Returns
    -------
    adata : annData
        The count data has been normalized and log-transformed with an offset of 1.
        The offset of 1 ensures that zero counts map to zeros. Keep this data in '.raw' part of the AnnData object
        as it will be used to visualize gene expression and perform statistical tests such as computing marker genes
        for clusters

    """

    # Keep the count data in a counts layer
    adata.layers["counts"] = adata.X.copy()

    # Normalize adata
    adata.X /= adata.obs['size_factors'].values[:, None]

    # log-transforms and updates adata
    # and log or Box-Cox transformation (lambda=0)
    # because count distribution follows a power-law and afterwards a normal distribution -> easier to apply stat-tests
    sc.pp.log1p(adata)

    # modify resulting matrix
    adata.X = np.asarray(adata.X)

    # Store the full data set in 'raw' as log-normalised data for statistical testing
    adata.raw = adata

    # save adata object
    sc.write('{}_QC_sizefactors.h5'.format(adata_path_filenames), adata=adata)

    return adata


def normalize_by_scanpy(adata, adata_path_filenames, exclude_highly_expressed=True, raw=False):
    """
    Normalize counts per spot (cell for scRNA-seq) with scanpy function sc.pp.normalize_total().
    If target sum is 1e6 than CPM normalisation is applied
    By excluding highly expressed genes, the normalisation of lower expressed genes are higher weighted [Weinreb17].

    Parameters
    ----------
    adata : annData
    adata_path_filenames : str
    exclude_highly_expressed : bool
    raw : bool

    Returns
    -------
    adata : annData
        The count data has been normalized and log-transformed with an offset of 1.
        The offset of 1 ensures that zero counts map to zeros. Keep this data in '.raw' part of the AnnData object
        as it will be used to visualize gene expression and perform statistical tests such as computing marker genes
        for clusters

    """
    # Keep the count data in a counts layer
    adata.layers["counts"] = adata.X.copy()

    # return dictionary if inplace is False otherwise updates adata
    x_norm = sc.pp.normalize_total(adata, inplace=False, exclude_highly_expressed=exclude_highly_expressed,
                                   target_sum=1e6)['X']
    adata.X = x_norm

    # log-transforms and updates adata
    sc.pp.log1p(adata)

    # modify resulting matrix
    adata.X = np.asarray(adata.X)

    # Store the full data set in 'raw' as log-normalised data for statistical testing
    adata.raw = adata

    # save adata object
    if raw:
        sc.write('{}_raw_QC_sizefactors.h5'.format(adata_path_filenames), adata=adata)
    else:
        sc.write('{}_QC_sizefactors.h5'.format(adata_path_filenames), adata=adata)

    return adata
