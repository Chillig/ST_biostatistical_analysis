from python_scripts.pre_processing.plot_functions import plots_preprocessing

import numpy as np
import scanpy as sc


def qc_ncounts_ngenes(adata, counts_threshold, save_folder, sample_name, lower_filter_counts=400,
                      upper_filter_counts=4000, upper_filter_genes=500, log_scale=None, raw=False):
    """Quality control - calculate QC covariates

    Parameters
    ----------
    adata : annData
        contains features and barcodes (in best: 1 barcode => 1 cell)
    counts_threshold : int
    save_folder : str
    sample_name : str
    lower_filter_counts : int
    upper_filter_counts : int
    upper_filter_genes : int
    log_scale : bool
    raw : bool

    Returns
    -------
    adata : annData
        contains now observables n_counts, log_counts, n_genes and mt_frac

    """
    # 1. cell QC
    # check for:
    #   number of molecule counts (UMIs)
    #   number of expressed genes
    #   fraction of spot with high mitochondrial percentage

    # sum up the number of UMI counts (barcodes) in count matrix .X along cols (genes)
    # number of counts per spot
    adata.obs['n_counts'] = adata.X.sum(1)
    # number of counts per spot in log
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    # number of genes per spot
    adata.obs['n_genes'] = (adata.X > 0).sum(1)

    mt_gene_mask = [gene.startswith('MT-') for gene in adata.var_names]
    adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1) / adata.obs['n_counts']

    # Threshold determination of UMI counts and genes
    qc_threshold_determination(adata=adata, counts_threshold=counts_threshold, save_folder=save_folder,
                               sample_name=sample_name, lower_filter_counts=lower_filter_counts,
                               upper_filter_counts=upper_filter_counts, upper_filter_genes=upper_filter_genes,
                               log_scale=log_scale, raw=raw)

    return adata


def qc_threshold_determination(adata, counts_threshold, save_folder, sample_name, lower_filter_counts=400,
                               upper_filter_counts=4000, upper_filter_genes=500, log_scale=None, raw=False):
    """Determine UMI count and gene threshold based on plots

    Parameters
    ----------
    adata : annData
    counts_threshold : int
    save_folder : str
    sample_name : str
    lower_filter_counts : int
    upper_filter_counts : int
    upper_filter_genes : int
    log_scale : bool
    raw : bool

    Returns
    -------

    """
    # Quality control - plot QC metrics
    plots_preprocessing.plot_qc_metrics(adata, save_folder=save_folder, sample_name=sample_name,
                                        counts_threshold=counts_threshold, raw=raw)

    # Thresholding decision: counts and genes
    plots_preprocessing.plot_distribution(adata, save_folder=save_folder, sample_name=sample_name,
                                          lower_filter_counts=lower_filter_counts,
                                          upper_filter_counts=upper_filter_counts,
                                          upper_filter_genes=upper_filter_genes,
                                          log_scale=log_scale, bins=[60, 25, 60, 60, 20], raw=raw)


def qc_filter(adata, min_counts_gene, max_counts_gene, mt_threshold, min_genes, min_shared_counts,
              max_gene_count_threshold=False, min_counts_spots=None, max_counts_spots=None):
    """Filter count matrix according to identified and given QC thresholds --> sort out outliers

    Parameters
    ----------
    adata : annData
    min_counts_gene : int
    max_counts_gene : int
    mt_threshold : float
    min_genes : int
    min_shared_counts : int
    max_gene_count_threshold : bool
        determine if in spatial transcriptomics max gene count threshold is useful
    min_counts_spots : int
    max_counts_spots : int

    Returns
    -------
    adata : annData
        filtered adata object

    """
    # Filter MT-fraction:
    print('Total number of spots: {:d}'.format(adata.n_obs))
    # Threshold for MT-fraction is 20% to 25%
    # ==> filter out spots with an overall high proportion of dying or highly stressed cells
    adata = adata[adata.obs['mt_frac'] < mt_threshold]
    print('Number of spots after MT filter: {:d}'.format(adata.n_obs))

    # Filter spots:
    # Minimum number of counts required for a spot to pass filtering.
    sc.pp.filter_cells(adata, min_counts=min_counts_spots)
    print('Number of spots after min count filter: {:d}'.format(adata.n_obs))

    # Maximum number of counts required for a spot to pass filtering.
    sc.pp.filter_cells(adata, max_counts=max_counts_spots)
    print('Number of spots after max count max_counts_spots: {:d}'.format(adata.n_obs))

    # Filter out spots which have a low number of genes expressed
    sc.pp.filter_cells(adata, min_genes=min_genes)
    print('Number of spots after gene filter: {:d}'.format(adata.n_obs))

    # Filter genes:
    print('Total number of genes: {:d}'.format(adata.n_vars))

    # Minimum number of cells expressed required for a gene to pass filtering
    sc.pp.filter_genes(adata, min_cells=min_shared_counts)
    print('Number of genes after spots filter: {:d}'.format(adata.n_vars))

    # Min. 20 UMI-counts - filters out 0 count genes by minimum number of counts required for a gene to pass filtering.
    sc.pp.filter_genes(adata, min_counts=min_counts_gene)

    # --> we have bulk resolution
    if max_gene_count_threshold:
        # Filters out genes by maximum number of counts required for a gene to pass filtering.
        sc.pp.filter_genes(adata, min_counts=max_counts_gene)

    return adata
