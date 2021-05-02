#! /usr/bin/python
"""Start pre-processing of data
    File name: main_preprocessing.py
    Author: Christina Hillig
    Date created: April/02/2020
    Date last modified: May/02/2021
    Python Version: 3.7
"""
from python_scripts.pre_processing import batch_correction, cell_cycle_storing, pc_determination, scaling_and_regression, \
    calculate_visualizations, quality_control, normalization, highly_variable_genes, init_variables, doublet_detection
from python_scripts.utils import loading_matrices, helper_tools as ht, sc_loading_matrices
from python_scripts.pre_processing.plot_functions import plots_preprocessing
import python_scripts.pre_processing.plot_functions.plot_imagespots as im_pp_plot

import scanpy as sc
import numpy as np
import os
from datetime import date


def load_dataset(configs):
    """Read out data sets, load images and spot locations in images

    Parameters
    ----------
    configs : configparser
        config file

    Returns
    -------
    raw_adata : annData
    filtered_adata : annData
    configs_path : str
    list_sample_ids : list

    """
    if configs['data']['data_type'] == 'Spatial Transcriptomics':
        raw_adata, filtered_adata, configs_path, list_sample_ids = loading_matrices.main(
            filename=configs['input_files']['rawdata_paths'],
            read_raw_matrix=configs.getboolean('preprocessing', 'read_raw_matrix'),
            spatial_concat=configs.getboolean('preprocessing', 'sample_concat'))
    else:
        raw_adata, filtered_adata, configs_path = sc_loading_matrices.main(
            configs['input_files']['rawdata_paths'],
            read_raw_matrix=configs.getboolean('preprocessing', 'read_raw_matrix'))
        list_sample_ids = filtered_adata.obs['sample'].values

    return raw_adata, filtered_adata, configs_path, list_sample_ids


def determine_heg(adata, save_folder):
    """Plot the highest expressed genes in a sample or whole data set
    Genes which are highly expressed can be, in some (disease) contexts, biologically relevant

    Parameters
    ----------
    adata : annData
    save_folder : str

    Returns
    -------

    """
    plots_preprocessing.plot_highest_expr_genes(adata=adata, save_folder=save_folder)


def plot_images(configs, adata, save_folder):
    """Visualise spots as grey dots on H&E slice

    Parameters
    ----------
    configs : configparser
    adata : annData
    save_folder : str
        path to result folder

    Returns
    -------

    """
    im_pp_plot.plot_greyspots_image(configs=configs, adata=adata, save_folder=save_folder, label="")


def sample_qc(adata, save_folder, sample_name, counts_threshold=60000, lower_filter_counts=400,
              upper_filter_counts=2500, upper_filter_genes=2000, log_scale=False, raw=False):
    """Check the quality of each sample and apply standard thresholds only for visualisation

    Parameters
    ----------
    adata : annData
    save_folder : str
        path to results folder
    sample_name : str
        name of sample or whole dataset
    counts_threshold : int
        threshold for UMI-counts
    lower_filter_counts : int
    upper_filter_counts : int
    upper_filter_genes : int
    log_scale : bool
        if plots should be shown in log-scale
    raw : bool
        if to load raw matrix from 10x Genomics

    Returns
    -------

    """
    adata_qc = quality_control.qc_ncounts_ngenes(
        adata, save_folder=save_folder, sample_name=sample_name,
        counts_threshold=counts_threshold, lower_filter_counts=lower_filter_counts,
        upper_filter_counts=upper_filter_counts, upper_filter_genes=upper_filter_genes,
        log_scale=log_scale, raw=raw)

    return adata_qc


def apply_qc_filter(adata, apply_mt_threshold):
    """Threshold determination of MT-fraction, UMI counts and genes per spot and per gene in all spots

    Parameters
    ----------
    adata : annData
    apply_mt_threshold : bool
        if a MT-fraction threshold shall be applied

    Returns
    -------

    """
    # Filter genes
    try:
        min_genes = int(input("Please provide the minimal No. genes threshold: "))
    except ValueError:
        print("Minimal No. genes threshold is now set to 30")
        min_genes = int(30)
    try:
        min_shared_counts = int(input("Please provide the minimal No. spots a gene has to be expressed; default 20: "))
    except ValueError:
        print("Minimal No. spots is now set to 20")
        min_shared_counts = int(20)
    try:
        min_umi_genes = int(input("Please provide the minimal No. UMI counts for a gene to pass filtering: "))
    except ValueError:
        print("Minimal No. UMI counts threshold is now set to 10")
        min_umi_genes = int(10)

    # Filter spots
    # min_counts=20, max_counts=60000, min_genes=20, min_cells=20)
    try:
        min_counts = int(input("Please provide the minimal count threshold: "))
    except ValueError:
        print("Minimal No. UMI counts threshold is now set to 5000")
        min_counts = int(5000)
    try:
        max_counts = int(input("Please provide the maximal count threshold: "))
    except ValueError:
        print("Minimal No. UMI counts threshold is now set to 35000")
        max_counts = int(35000)
    # ------------------------------------------------------------------------------------------------------------ #
    if apply_mt_threshold:
        """ATTENTION: Regressing out biological covariates is generally done to isolate particular processes in 
                        the data in which you are interested in, while losing global structure in the data.
        MT gene expression is also a biological covariate (as well as a technical indicator of cell stress) """
        try:
            mt_threshold = float(input("Please provide MT-fraction threshold; recom. 0.2 - 0.25 (default 0.25): "))
        except ValueError:
            print("MT-fraction threshold is now set to 25%")
            mt_threshold = float(0.25)
    else:
        # maximum value = 100% MT-fraction
        print("MT-fraction threshold is now set to 100%")
        mt_threshold = float(1.0)
    # ------------------------------------------------------------------------------------------------------------ #

    gene_count_threshold = input("Apply max UMI-threshold for a gene to pass filtering (y/n)?: ")
    if gene_count_threshold == 'y':
        # TODO investigate if it makes sense to filter for max UMI counts...
        #  (eg outliers should be removed by normalisation)
        try:
            max_umi_genes = int(input("Please provide the maximal No. UMI counts for a gene to pass filtering: "))
        except ValueError:
            print("Maximal No. UMI counts threshold is now set to 5000")
            max_umi_genes = int(5000)

        cutted_adata = quality_control.qc_filter(adata,
                                                 min_counts_spots=min_counts, max_counts_spots=max_counts,
                                                 min_genes=min_genes,
                                                 min_counts_gene=min_umi_genes, max_counts_gene=max_umi_genes,
                                                 min_shared_counts=min_shared_counts,
                                                 mt_threshold=mt_threshold,
                                                 max_gene_count_threshold=True)

    else:
        max_umi_genes = 0
        cutted_adata = quality_control.qc_filter(adata, min_counts_spots=min_counts, max_counts_spots=max_counts,
                                                 min_genes=min_genes, min_counts_gene=min_umi_genes,
                                                 max_counts_gene=max_umi_genes,
                                                 min_shared_counts=min_shared_counts,
                                                 mt_threshold=mt_threshold, max_gene_count_threshold=False)

    return \
        cutted_adata, min_genes, min_shared_counts, mt_threshold, min_counts, max_counts, min_umi_genes, max_umi_genes


def apply_normalisation(adata, save_folder, adata_path_filenames, norm_type="scanpy", exclude_highly_expressed=True,
                        raw=False):
    """Normalize the raw counts to account for differences in sequencing depth per cell for each sample

    Parameters
    ----------
    adata : annData
    save_folder : str
    adata_path_filenames : str
    norm_type : str
    exclude_highly_expressed : bool
    raw : bool

    Returns
    -------

    """
    if norm_type == "scran":
        # Perform a clustering for scran normalization in clusters but takes much longer compared to scanpy norm method
        adata = normalization.first_est_sizefactors(adata=adata, save_folder=save_folder, raw=raw)
        # # Normalize data by sizefactor
        norm_adata = normalization.normalize_by_sizefactor(adata=adata, adata_path_filenames=adata_path_filenames)
    else:
        norm_adata = normalization.normalize_by_scanpy(adata=adata, adata_path_filenames=adata_path_filenames,
                                                       exclude_highly_expressed=exclude_highly_expressed, raw=raw)

    return norm_adata


def add_metadata(adata):
    """Add metaData as observable to annData object
        -> data specific

    Parameters
    ----------
    adata : annData

    Returns
    -------

    """
    # get labels
    tissue_cell_labels, disease_labels, lesion_labels = ht.get_tissue_annot(adata)
    # get batches assigned to each patient (currently we have 2 or 4 samples per patient)
    df_batches = ht.map_sample_batch_list(adata=adata,
                                          num_samples_patient=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2])
    # assign donor to each spot
    donors = []
    for project in df_batches:
        for ind, s_c in enumerate(df_batches[project]['sample']):
            no_spots = adata[adata.obs['sample'] == s_c].shape[0]
            donors.extend(np.ones(no_spots) * df_batches[project]['batch'][ind])
    # need information about how often samples were found
    adata.obs['patient'] = donors
    adata.obs['patient'] = adata.obs['patient'].astype('int').astype('category')

    # assign Diagnostics
    spot_disease = []
    for spot_label in disease_labels:
        spots = adata[adata.obs[spot_label] == 1]
        spot_disease.extend([spot_label] * spots.shape[0])
    adata.obs['disease'] = spot_disease
    adata.obs['disease'] = adata.obs['disease'].astype('string').astype('category')

    # assign Biopsy type
    spot_biopsy = []
    for spot_label in lesion_labels:
        spots = adata[adata.obs[spot_label] == 1]
        spot_biopsy.extend([spot_label] * spots.shape[0])
    adata.obs['biopsy_type'] = spot_biopsy
    adata.obs['biopsy_type'] = adata.obs['biopsy_type'].astype('string').astype('category')

    return adata, tissue_cell_labels, disease_labels, lesion_labels


def main(configs, adata, save_folder):
    """Control which pre-processing steps shall be applied before downstream analysis, DGE analysis
    and pathway enrichment analysis

    Parameters
    ----------
    configs : configparser
        contains all parameters -> to add: thresholds and cut parameters
    adata : annData
    save_folder : str

    Returns
    -------
    norm_pp_adata : annData
    adata_filename : str

    """

    print("\n-------- Overview of data sets --------")
    # 1.1 Add meta data like which samples belong to which donor (since 04.10.2020)
    adata, tissue_cell_labels, disease_labels, lesion_labels = add_metadata(adata)
    # 1.2 Remove spots having no tissue/cell labels (since 06.10.2020)
    adata = adata[np.where(adata.obs[tissue_cell_labels].to_numpy().any(axis=1))[0]]

    # print info about sample 1
    sample_name = adata.obs['sample'].values[1]
    print("\nSample {} ".format(sample_name))
    print("Shape of filtered data set: ", adata.shape)
    print("Tissue associated No. of spots: ", adata.shape[0])
    print("Tissue associated No. of genes: ", adata.shape[1])
    print("Observables contained in data sets sorted by barcodes: ", adata.obs_keys())
    print("Variables contained in data sets sorted by gene names: ", adata.var_keys())

    # plot spots on top of images (only for the first sample)
    plot_images(configs=configs, adata=adata, save_folder=save_folder)

    # 2. Pre-processing and visualization
    # apply the following steps 2.1 - 2.6 individually on each adata object
    print("\n-------- Start Pre-processing and Visualization --------")

    # 2.0
    # show 20thst highest expressed genes (HEG) in data set and per sample
    determine_heg(adata=adata, save_folder=save_folder)

    print("\n         Quality Control\n")
    # 2.1 QC (Quality control) of data - calculate QC covariates
    # 2.1.1 Cell QC
    # TODO Determine counts_threshold via Mean absolute deviation (MAD); find outliers :)
    adata_qc = sample_qc(adata=adata, save_folder=save_folder, counts_threshold=60000, lower_filter_counts=400,
                         upper_filter_counts=2500, upper_filter_genes=2000, log_scale=False,
                         raw=configs.getboolean("preprocessing", "read_raw_matrix"), sample_name=sample_name)

    # 2.1.2 Threshold determination of UMI counts and genes
    cutted_adata, min_genes, min_shared_counts, mt_threshold, min_counts, max_counts, min_umi_genes, max_umi_genes = \
        apply_qc_filter(adata=adata_qc,
                        apply_mt_threshold=configs.getboolean("preprocessing", 'apply_mt_threshold'))

    if configs['preprocessing']['filter_doublets']:
        # 2.1.3 Filter out multiplets
        cutted_adata = doublet_detection.scrublet_algorithm(cutted_adata, save_folder=save_folder)

    # save QC adata object
    sc.write('{}_QC.h5'.format(configs["data"]['output_path']), cutted_adata)

    # 2.2 Normalization
    print("\n         Normalization\n")
    norm_adata = apply_normalisation(adata=cutted_adata, save_folder=save_folder,
                                     norm_type=configs['preprocessing']['normalisation_function'],
                                     exclude_highly_expressed=configs.getboolean("preprocessing",
                                                                                 "exclude_highly_expressed"),
                                     raw=configs_file.getboolean("preprocessing", "read_raw_matrix"),
                                     adata_path_filenames=configs["data"]['output_path'])
    # TODO plot normalised count data distribution

    # -------------------------------------------------- Optional ---------------------------------------------------- #
    # 2.2.1 Scale data
    if configs.getboolean("preprocessing", "apply_scaling"):
        norm_adata = scaling_and_regression.scaling(adata=norm_adata)

    # 2.2.2 Regress out uninteresting variation
    """ATTENTION: Regressing out biological covariates is generally done to isolate particular processes in 
                  the data that you are interested in, while losing global structure in the data.
    Cell cycle stage can be a major determinant in the difference between two cell types (e.g. stem cells and 
    proliferating cells like transit amplifying cells). Removing this effect, hides the distinction"""
    if configs.getboolean("preprocessing", "apply_remove_cc_effect"):
        norm_adata = scaling_and_regression.apply_regression_variables(adata=norm_adata)
    # ---------------------------------------------------------------------------------------------------------------- #

    print("\n         Visualisation\n")
    # 2.3.1 Highly Variable Genes (HVG) for feature selection
    # HVG: highly expressed in some cells and lowly expressed in others
    norm_adata = highly_variable_genes.find_highly_variable_genes(
        norm_adata, type_dataset="uncorrected", save_folder=save_folder, num_top_genes=4000,
        raw=configs.getboolean("preprocessing", "read_raw_matrix"))

    # 2.3.2 Determine No. PCs
    pc_determination.pcs_combs(norm_adata, save_folder, use_highly_variable=True, copy=False, return_info=False,
                               raw=configs.getboolean("preprocessing", "read_raw_matrix"), type_dataset="uncorrected")

    # 2.3.3 Visualization
    try:
        n_comps = int(input("Please provide the No. principal components (default 50): "))
    except ValueError:
        n_comps = int(50)
    norm_adata = calculate_visualizations.calc_visualization(
        norm_adata, save_folder=save_folder, raw=configs.getboolean("preprocessing", "read_raw_matrix"),
        n_comps=n_comps, batch_key="uncorrected")

    if configs.getboolean("preprocessing", "get_cc_effect"):
        print("\n         Cell Cycle\n")
        # 2.4 Cell cycle scoring
        # todo find another .csv file with cell cycle phases
        norm_adata = cell_cycle_storing.score_cell_cycle(
            cc_genes_file=configs['input_files']['cell_cycle'], adata=norm_adata, save_folder=save_folder,
            raw=configs.getboolean("preprocessing", "read_raw_matrix"))

    # 2.5 Apply Batch correction if samples are from same (or different) data set but splitted into batches
    # Dr. Maren Buettner:
    # "During the QC step, we observed differences across samples for instance, in the library size per dataset.
    # Such differences may contribute to the batch effect."
    if configs.getboolean("preprocessing", "sample_concat"):
        print("\n         Batch Correction\n")
        norm_bc_adata = batch_correction.apply_batch_correction(norm_adata, save_folder=save_folder, n_comps=n_comps)

        # 2.5.1 Run find highly variable genes again on integrated dataset
        # HVG: highly expressed in some cells and lowly expressed in others
        norm_pp_adata = highly_variable_genes.find_highly_variable_genes(
            norm_bc_adata, type_dataset="batch_corrected", save_folder=save_folder, num_top_genes=4000,
            raw=configs.getboolean("preprocessing", "read_raw_matrix"))
        # Actually already calculated in Batch correction functions..
        # 2.5.2 Determine No. PCs
        pc_determination.pcs_combs(norm_pp_adata, save_folder, type_dataset="batch_corrected",
                                   raw=configs.getboolean("preprocessing", "read_raw_matrix"))

        # 2.5.3 Visualisation
        n_comps = int(input("Please provide the No. principal components (default 50): "))
        norm_pp_adata = calculate_visualizations.calc_visualization(
            norm_pp_adata, save_folder=save_folder, raw=configs.getboolean("preprocessing", "read_raw_matrix"),
            n_comps=n_comps, batch_key="batch_corrected")

    else:
        norm_pp_adata = norm_adata.copy()

    sc.write('{}_QC_BC.h5'.format(configs["data"]['output_path']), norm_pp_adata)
    plots_preprocessing.plot_visualization_results(
        adata=norm_pp_adata, save_folder=save_folder, batch_key="batch_corrected",
        raw=configs.getboolean("preprocessing", "read_raw_matrix"))

    print("-------- Finished: Pre-processing and Visualization --------")

    print("Start storing pre-processed AnnData object")
    # 2.7 save pre-processed annData object
    # # transform float e.g. 0.25 -> 0_25
    mt_cut_splitted = str(mt_threshold).split(".")
    mt_cut = mt_cut_splitted[0] + str("_") + mt_cut_splitted[1]

    if max_umi_genes == 0:
        # save pre-processed annData object
        filter_name = 'minumi_{}_maxumi_{}_mg_{}_msc_{}_mt_{}_minumig{}'.format(min_counts, max_counts, min_genes,
                                                                                min_shared_counts, mt_cut,
                                                                                min_umi_genes)
    else:
        # save pre-processed annData object
        filter_name = 'minumi_{}_maxumi_{}_mg_{}_msc_{}_mt_{}_minumig{}_manumig{}'.format(min_counts, max_counts,
                                                                                          min_genes, min_shared_counts,
                                                                                          mt_cut, min_umi_genes,
                                                                                          max_umi_genes)

    adata_filename = '{}_{}_pp.h5'.format(configs["data"]['output_path'], filter_name)
    sc.write(adata_filename, norm_pp_adata)

    return norm_pp_adata, adata_filename


if __name__ == '__main__':
    output_path = os.path.join("..", "..", "output", str(date.today()))
    adata_savepath = init_variables.init_vars()
    configs_file = ht.load_config(config_path=adata_savepath)

    print("\nCapture area in Visium slide contains a grid of 4,992 capture spots")
    # 1. Load data
    print("#   --  >Load data and information<  --   #")
    _, unpp_filtered_adata, _, _ = load_dataset(configs=configs_file)

    print("-------- Finished: Read out values --------")
    pp_adata, filename_adata = main(configs=configs_file, adata=unpp_filtered_adata, save_folder=output_path)
