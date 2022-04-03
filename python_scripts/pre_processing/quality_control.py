from python_scripts.pre_processing.plot_functions import plots_preprocessing

import numpy as np
import scanpy as sc
import anndata
import configparser


def calculate_ncounts_ngenes(adata, key=''):
    """

    Parameters
    ----------
    adata : annData
    key : str

    Returns
    -------
    adata : annData
        contains now observables n_counts, log_counts, n_genes and mt_frac

    """
    # sum up the number of UMI counts (barcodes) in count matrix .X along cols (genes)
    # number of counts per spot
    adata.obs['n_counts{}'.format(key)] = adata.X.sum(1)
    # number of counts per spot in log
    adata.obs['log_counts{}'.format(key)] = np.log(adata.obs['n_counts{}'.format(key)])
    # number of genes per spot
    adata.obs['n_genes{}'.format(key)] = (adata.X > 0).sum(1)

    mt_gene_mask = [gene.startswith('MT-') for gene in adata.var_names]
    adata.obs['mt_frac{}'.format(key)] = adata.X[:, mt_gene_mask].sum(1) / adata.obs['n_counts{}'.format(key)]

    # add ribosomal genes
    # mitochondrial genes
    adata.var['mt{}'.format(key)] = adata.var_names.str.startswith('MT-')
    # ribosomal genes
    adata.var['ribo{}'.format(key)] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var['hb{}'.format(key)] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt{}'.format(key), 'ribo{}'.format(key), 'hb{}'.format(key)],
                               percent_top=None, log1p=True, inplace=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt{}'.format(key), 'ribo{}'.format(key), 'hb{}'.format(key)],
                               percent_top=None, log1p=False, inplace=True)

    print("Overview of n_counts: ", adata.obs['n_counts{}'.format(key)].describe())
    print("Overview of n_genes: ", adata.obs['n_genes{}'.format(key)].describe())

    return adata


def qc_ncounts_ngenes(adata: anndata, save_folder: str, configs: configparser,
                      counts_threshold: int = 60000, lower_filter_counts: int = 400, upper_filter_counts: int = 2500,
                      upper_filter_genes: int = 2000, log_scale: bool = False):
    """Quality control - calculate QC covariates

    Parameters
    ----------
    adata : annData
        contains features and barcodes (in best: 1 barcode => 1 cell)
    configs : ConfigParser
    counts_threshold : int
        threshold for UMI-counts
    save_folder : str
        path to results folder
    lower_filter_counts : int
    upper_filter_counts : int
    upper_filter_genes : int
    log_scale : bool
        if plots should be shown in log-scale

    Returns
    -------
    adata : annData
        contains now observables n_counts, log_counts, n_genes and mt_frac

    """
    adata = calculate_ncounts_ngenes(adata)

    # Threshold determination of UMI counts and genes:
    # -> Determine UMI count and gene threshold based on plots
    # Quality control - plot QC metrics
    # TODO reactivate
    # plots_preprocessing.plot_qc_metrics(
    #     adata, save_folder=save_folder, project_name=configs['project']['project_name'],
    #     raw=configs.getboolean("preprocessing", "read_raw_matrix"))
    #
    # # Thresholding decision: counts and genes
    # plots_preprocessing.plot_distribution(adata, save_folder=save_folder,
    #                                       project_name=configs['project']['project_name'],
    #                                       lower_filter_counts=lower_filter_counts,
    #                                       upper_filter_counts=upper_filter_counts,
    #                                       upper_filter_genes=upper_filter_genes,
    #                                       log_scale=log_scale, bins=[60, 25, 60, 60, 20],
    #                                       raw=configs.getboolean("preprocessing", "read_raw_matrix"))
    #
    # plots_preprocessing.plot_project_counts_dist(
    #     adata=adata, save_folder=save_folder, raw=configs.getboolean("preprocessing", "read_raw_matrix"))

    return adata


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

    # --> we have bulk resolution
    if max_gene_count_threshold:
        # Filters out genes by maximum number of counts required for a gene to pass filtering
        # -> removes potential doublets
        sc.pp.filter_cells(adata, max_counts=max_counts_gene)

    # Filter MT-fraction:
    print('Total number of spots: {:d}'.format(adata.n_obs))
    # Threshold for MT-fraction is 20% to 25%
    # ==> filter out spots with an overall high proportion of dying or highly stressed cells
    adata = adata[adata.obs['mt_frac'] < mt_threshold]
    print('Number of spots after MT filter: {:d}'.format(adata.n_obs))

    # Filter genes:
    print('Total number of genes: {:d}'.format(adata.n_vars))

    # Minimum number of cells expressed required for a gene to pass filtering
    sc.pp.filter_genes(adata, min_cells=min_shared_counts)
    print('Number of genes after spots filter: {:d}'.format(adata.n_vars))

    # Min. 20 UMI-counts - filters out 0 count genes by minimum number of counts required for a gene to pass filtering.
    sc.pp.filter_genes(adata, min_counts=min_counts_gene)

    return adata
