#! /usr/bin/python
"""Start pre-processing of data
    File name: main_preprocessing.py
    Author: Christina Hillig
    Date created: April/02/2020
    Date last modified: May/02/2021
    Python Version: 3.7
"""
from python_scripts.pre_processing import batch_correction, cell_cycle_storing, pc_determination, \
    scaling_and_regression, calculate_visualizations, quality_control, normalization, highly_variable_genes, \
    init_variables, doublet_detection
from python_scripts.utils import loading_matrices, helper_tools as ht, sc_loading_matrices, add_observables
from python_scripts.pre_processing.plot_functions import plots_preprocessing
import python_scripts.pre_processing.plot_functions.plot_imagespots as im_pp_plot

import scanpy as sc
import numpy as np
import pandas as pd
import os
from datetime import date
import configparser
import glob


def load_dataset(configs):
    """Read out data sets, load images and spot locations in images

    Parameters
    ----------
    configs : ConfigParser
        config file

    Returns
    -------
    raw_adata : annData
    filtered_adata : annData
    configs_path : str
    list_sample_ids : list

    """
    if configs['data']['data_type'] == 'Spatial Transcriptomics':
        print("\nCapture area in Visium slide contains a grid of 4,992 capture spots")
        raw_adata, filtered_adata, configs_path, list_sample_ids = loading_matrices.main(
            filename=configs['input_files']['rawdata_paths'],
            read_raw_matrix=configs.getboolean('preprocessing', 'read_raw_matrix'),
            spatial_concat=configs.getboolean('preprocessing', 'sample_concat'))
    else:
        raw_adata, filtered_adata, configs_path = sc_loading_matrices.main(
            filename=configs['input_files']['rawdata_paths'],
            read_raw_matrix=configs.getboolean('preprocessing', 'read_raw_matrix'))
        list_sample_ids = filtered_adata.obs['sample'].values

    return raw_adata, filtered_adata, configs_path, list_sample_ids


def determine_heg(adata, save_folder):
    """Plot the highest expressed genes in a Disease or whole data set
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
    # ------------------------------------------------------------------------------------------------------------ #
    if apply_mt_threshold:
        print("Filter MT-fraction")
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

    # Filter genes
    print("Filter genes")
    try:
        min_genes = int(input("Please provide the minimal No. genes threshold; default 30: "))
    except ValueError:
        print("Minimal No. genes threshold is now set to 30")
        min_genes = int(30)
    try:
        min_shared_counts = int(input("Please provide the minimal No. spots a gene has to be expressed; default 20: "))
    except ValueError:
        print("Minimal No. spots is now set to 20")
        min_shared_counts = int(20)
    try:
        min_umi_genes = int(
            input("Please provide the minimal No. UMI counts for a gene to pass filtering; default 1: "))
    except ValueError:
        print("Minimal No. UMI counts threshold is now set to 1")
        min_umi_genes = int(1)

    # Filter spots
    print("Filter cells")
    # min_counts=20, max_counts=60000, min_genes=20, min_cells=20)
    try:
        min_counts = int(input("Please provide the minimal count threshold; default 500: "))
    except ValueError:
        print("Minimal No. UMI counts threshold is now set to 500")
        min_counts = int(500)
    try:
        max_counts = int(input("Please provide the maximal count threshold; default 35000: "))
    except ValueError:
        print("Maximal No. UMI counts threshold is now set to 35000")
        max_counts = int(35000)

    gene_count_threshold = input("Apply max UMI-threshold for a gene to pass filtering (y/n)?: ")
    if gene_count_threshold == 'y':
        # TODO investigate if it makes sense to filter for max UMI counts...
        #  (eg outliers should be removed by normalisation)
        try:
            max_umi_genes = int(
                input("Please provide the maximal No. UMI counts for a gene to pass filtering; default 5000: "))
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


def apply_normalisation(adata, save_folder, configs):
    """Normalize the raw counts to account for differences in sequencing depth per cell for each specimen

    Parameters
    ----------
    adata : annData
    save_folder : str
    configs : ConfigParser

    Returns
    -------

    """
    if configs['preprocessing']['normalisation_function'] == "scran":
        # Perform a clustering for scran normalization in clusters but takes much longer compared to scanpy norm method
        adata = normalization.first_est_sizefactors(adata=adata, save_folder=save_folder, configs=configs)
        # # Normalize data by sizefactor
        norm_adata = normalization.normalize_by_sizefactor(adata=adata, configs=configs)
    else:
        norm_adata = normalization.normalize_by_scanpy(adata=adata, configs=configs)

    return norm_adata


def main(configs, adata, save_folder):
    """Control which pre-processing steps shall be applied before downstream analysis, DGE analysis
    and pathway enrichment analysis

    Parameters
    ----------
    configs : ConfigParser
        contains all parameters -> to add: thresholds and cut parameters
    adata : annData
    save_folder : str

    Returns
    -------
    norm_pp_adata : annData
    adata_filename : str

    """

    print("\n-------- Overview of data sets --------")
    if configs['data']['data_type'] == 'Spatial Transcriptomics':
        dataset_type = "st"
    else:
        dataset_type = "sc"

    adata = adata[adata.obs['DISEASE'].str.contains('Pso|AD|PRP')].copy()

    # print info about sample 1
    print("\nDataset {} ".format(dataset_type))
    print("Shape of filtered data set: ", adata.shape)
    print("Tissue associated No. of spots: ", adata.shape[0])
    print("Tissue associated No. of genes: ", adata.shape[1])
    print("Observables contained in data sets sorted by barcodes: ", adata.obs_keys())
    print("Variables contained in data sets sorted by gene names: ", adata.var_keys())

    # if configs['data']['data_type'] == 'Spatial Transcriptomics':
    #     # plot spots on top of images (only for the first sample)
    #     plot_images(configs=configs, adata=adata, save_folder=save_folder)
    # sc.pl.spatial(adata_prp[adata_prp.obs['Sample_CaptureArea'] == 'P21093_21L008977'],
    #               color='patient', size=1, library_id='P21093_21L008977')

    # 2. Pre-processing and visualization
    # apply the following steps 2.1 - 2.6 individually on each adata object
    print("\n-------- Start Pre-processing and Visualization --------")

    # 2.0.1
    # show 20thst highest expressed genes (HEG) in data set and per DISEASE
    determine_heg(adata=adata, save_folder=save_folder)

    # 2.0.2
    if configs.getboolean('preprocessing', 'filter_doublets'):
        # Check doublets prior to QC
        adata = doublet_detection.carlos_woublet(adata=adata, save_folder=save_folder, key='pre_QC')

    print("\n         Quality Control\n")
    # 2.1 QC (Quality control) of data - calculate QC covariates
    adata_filename = '{}_QC.h5'.format(dataset_type)
    # 2.1.1 Cell QC
    if not os.path.exists(os.path.join(configs["data"]['output_path'], adata_filename)):
        # TODO Determine counts_threshold via Mean absolute deviation (MAD); find outliers :)
        adata_qc = quality_control.qc_ncounts_ngenes(
            adata=adata, configs=configs, save_folder=save_folder,
            counts_threshold=80000, lower_filter_counts=2000,
            upper_filter_counts=2500, upper_filter_genes=700, log_scale=False)

        # 2.1.2 Threshold determination of UMI counts and genes
        cutted_adata, min_genes, min_shared_counts, mt_threshold, min_counts, max_counts, \
        min_umi_genes, max_umi_genes = apply_qc_filter(
            adata=adata_qc, apply_mt_threshold=configs.getboolean("preprocessing", 'apply_mt_threshold'))

        # add cut-offs to config file
        config = configparser.ConfigParser()
        if not config.has_section("Set_parameters"):
            config.add_section("Set_parameters")

            # add section using append
            section = 'a'
        else:
            # overwrite parameters using wb
            section = 'wb'

        config.set("Set_parameters", "min_genes", str(min_genes))
        config.set("Set_parameters", "min_shared_counts", str(min_shared_counts))
        config.set("Set_parameters", "mt_threshold", str(mt_threshold))
        config.set("Set_parameters", "min_counts", str(min_counts))
        config.set("Set_parameters", "max_counts", str(max_counts))
        config.set("Set_parameters", "min_umi_genes", str(min_umi_genes))
        if max_umi_genes == 0:
            max_umi_genes = np.Inf
        config.set("Set_parameters", "max_umi_genes", str(max_umi_genes))

        # Write parameters to config file
        configfile_list = glob.glob(os.path.join(adata_savepath, "*.ini"))
        with open(configfile_list[0], section) as configfile:
            config.write(configfile)

        # recalculate n_counts and n_genes
        cutted_adata = quality_control.calculate_ncounts_ngenes(cutted_adata, key='_qc')

        sc.write(os.path.join(configs["data"]['output_path'], adata_filename), cutted_adata)

        cutted_adata.obs['tissue_layer'].to_csv(os.path.join(configs["data"]['output_path'], 'Tissue_layers.csv'))
        cutted_adata.obs['spot_type'].to_csv(os.path.join(configs["data"]['output_path'], 'Spot_types.csv'))
    else:
        del adata
        cutted_adata = sc.read(os.path.join(configs["data"]['output_path'], adata_filename))

    if configs.getboolean('preprocessing', 'filter_doublets'):
        adata_filename = '{}_QC_DD.h5'.format(dataset_type)
        if not os.path.exists(os.path.join(configs["data"]['output_path'], adata_filename)):
            # 2.1.3 Filter out multiplets --
            cutted_adata = doublet_detection.carlos_woublet(adata=cutted_adata, save_folder=save_folder, key='post_QC')
            # cutted_adata = doublet_detection.scrublet_algorithm(
            #     cutted_adata, sample_name=sample_name, save_folder=save_folder)

            # save QC and adata object
            sc.write(os.path.join(configs["data"]['output_path'], adata_filename), cutted_adata)

    # 2.2 Normalization
    print("\n         Normalization\n")
    # TODO apply normalisation on skin layers and include both spot type and skin layers in design
    adata_filename = '{}_QC_normed.h5'.format(dataset_type)
    if not os.path.exists(os.path.join(configs["data"]['output_path'], adata_filename)):
        norm_adata = apply_normalisation(adata=cutted_adata, save_folder=save_folder, configs=configs)
        # save QC and normed adata object
        sc.write(os.path.join(configs["data"]['output_path'], adata_filename), norm_adata)
        # TODO plot normalised count data distribution
    else:
        del cutted_adata
        norm_adata = sc.read(os.path.join(configs["data"]['output_path'], adata_filename))

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

    # 2.5 Apply Batch correction if specimen are from same (or different) data set but splitted into batches
    # Dr. Maren Buettner:
    # "During the QC step, we observed differences across samples for instance, in the library size per dataset.
    # Such differences may contribute to the batch effect."
    if configs.getboolean("preprocessing", "apply_batch_correction"):
        print("\n         Batch Correction\n")
        norm_bc_adata = batch_correction.apply_batch_correction(
            norm_adata, save_folder=save_folder, n_comps=n_comps, batch_key='specimen',
            bc_tool=configs["preprocessing"]['batch_correction_tool'],
            possible_batch_effects=['sample', 'specimen', 'project', 'phase', 'patient', 'disease', 'biopsy_type',
                                    'capture_area', 'object_slide'])
        del norm_adata

        # 2.5.1 Run find highly variable genes again on integrated dataset
        # HVG: highly expressed in some cells and lowly expressed in others
        norm_pp_adata = highly_variable_genes.find_highly_variable_genes(
            norm_bc_adata, type_dataset="batch_corrected", save_folder=save_folder, num_top_genes=4000,
            raw=configs.getboolean("preprocessing", "read_raw_matrix"))
        del norm_bc_adata
        # Actually already calculated in Batch correction functions..
        # 2.5.2 Determine No. PCs
        pc_determination.pcs_combs(norm_pp_adata, save_folder, type_dataset="batch_corrected",
                                   use_highly_variable=True, copy=False, return_info=False,
                                   raw=configs.getboolean("preprocessing", "read_raw_matrix"))

        # 2.5.3 Visualisation
        n_comps = int(input("Please provide the No. principal components (default 50): "))
        norm_pp_adata = calculate_visualizations.calc_visualization(
            norm_pp_adata, save_folder=save_folder, raw=configs.getboolean("preprocessing", "read_raw_matrix"),
            n_comps=n_comps, batch_key="batch_corrected")

        adata_filename = '{}_QC_normed_BC.h5'.format(dataset_type, )
        sc.write(os.path.join(configs["data"]['output_path'], adata_filename), norm_pp_adata)

    else:
        norm_pp_adata = norm_adata.copy()

    plots_preprocessing.plot_visualization_results(
        adata=norm_pp_adata, save_folder=save_folder, batch_key="batch_corrected",
        raw=configs.getboolean("preprocessing", "read_raw_matrix"))

    print("-------- Finished: Pre-processing and Visualization --------")

    print("Start storing pre-processed AnnData object")
    # 2.7 save pre-processed annData object
    # # transform float e.g. 0.25 -> 0_25
    mt_cut_splitted = str(configs["data"]['mt_threshold']).split(".")
    mt_cut = mt_cut_splitted[0] + str("_") + mt_cut_splitted[1]

    if configs["data"]['max_umi_genes'] == np.inf:
        # save pre-processed annData object
        filter_name = '{}_minumi_{}_maxumi_{}_mg_{}_msc_{}_mt_{}_minumig_{}'.format(
            dataset_type, configs["data"]['min_counts'], configs["data"]['max_counts'], configs["data"]['min_genes'],
            configs["data"]['min_shared_counts'], mt_cut, configs["data"]['min_umi_genes'])
    else:
        # save pre-processed annData object
        filter_name = '{}_minumi_{}_maxumi_{}_mg_{}_msc_{}_mt_{}_minumig{}_maxumig_{}'.format(
            dataset_type, configs["data"]['min_counts'], configs["data"]['max_counts'], configs["data"]['min_genes'],
            configs["data"]['min_shared_counts'], mt_cut, configs["data"]['min_umi_genes'],
            configs["data"]['max_umi_genes'])

    adata_filename = '{}_pp.h5'.format(filter_name)
    sc.write(os.path.join(configs["data"]['output_path'], adata_filename), norm_pp_adata)

    return norm_pp_adata, adata_filename


if __name__ == '__main__':
    # 1. Visualise MT-fraction on H&E images
    # 2. Identify spatially variable genes
    output_path = os.path.join("..", "..", "output", str(date.today()))
    os.makedirs(output_path, exist_ok=True)
    adata_savepath = init_variables.init_vars()
    configs_file = ht.load_config(config_path=adata_savepath)

    # 1. Load data
    print("#   --  >Load data and information<  --   #")
    unppadata_filename = '{}_unpp.h5'.format(configs_file['data']['data_type'])
    if not os.path.isfile(os.path.join(adata_savepath, unppadata_filename)):
        _, unpp_filtered_adata, _, _ = load_dataset(configs=configs_file)  # 93373 x 20613
        sc.write(os.path.join(adata_savepath, unppadata_filename), unpp_filtered_adata)
    else:
        unpp_filtered_adata = sc.read(os.path.join(adata_savepath, unppadata_filename))

    unppadata_cleaned_filename = '{}_unpp_cleaned.h5'.format(configs_file['data']['data_type'])
    if not os.path.isfile(os.path.join(adata_savepath, unppadata_cleaned_filename)):
        # Add patient info for new samples
        # 1. Remove all samples with no SAMPLE annotation (0)
        m_unkown_samples = (unpp_filtered_adata.obs['SAMPLE'] == 0) | (unpp_filtered_adata.obs['SAMPLE'] == '0')
        unpp_filtered_adata = unpp_filtered_adata[~m_unkown_samples].copy()  # -> removed 9746 spots -> 83627 left

        # 2. Add sample ids to sample names
        unpp_filtered_adata.obs['SAMPLE'] = unpp_filtered_adata.obs['SAMPLE'].astype(str)
        m_samplenames = ~pd.to_numeric(unpp_filtered_adata.obs['SAMPLE'], errors='coerce').isna()
        unpp_filtered_adata.obs.loc[m_samplenames, 'SAMPLE'] = unpp_filtered_adata.obs.loc[
            m_samplenames, 'SAMPLE'].replace(r'\.0$', '', regex=True)

        # # 1.1 Add meta data like which samples belong to which donor (since 04.10.2020)
        unpp_filtered_adata, _lesion_labels = add_observables.add_lesion_metadata(unpp_filtered_adata)

        # Save annotations for normalisation in .csv file
        # Use annotations from pathologist instead of clusters
        # Correct for tissue layers and large tissue structures:
        unpp_filtered_adata = add_observables.add_tissuelayers_obs(adata=unpp_filtered_adata)
        # Add spottype == large tissue structures:
        # SEBACEOUS GLAND, SWEAT GLAND, MUSCLE, HAIRFOLLICLE, VESSEL,
        # DERMIS, upper EPIDERMIS, middle and basal EPIDERMIS, JUNCTION
        unpp_filtered_adata = add_observables.add_spottypes_obs(adata=unpp_filtered_adata)
        # 1.2 Remove spots having no tissue label: 81331 × 20613  neu: 83627
        # TODO with new annotations: 78902
        unpp_filtered_adata = unpp_filtered_adata[~(unpp_filtered_adata.obs['tissue_layer'] == 'Unknown')].copy()

        # adjust for renamed column disease to DISEASE
        # 1. add DISEASE to disease column
        unpp_filtered_adata.obs['DISEASE'] = unpp_filtered_adata.obs['DISEASE'].astype('str')
        # 2. rename PSO to Pso
        unpp_filtered_adata.obs = unpp_filtered_adata.obs.replace(
            to_replace={'DISEASE': 'PSO'}, value={'DISEASE': 'Pso'}, regex=True)
        unpp_filtered_adata.obs['DISEASE'] = unpp_filtered_adata.obs['DISEASE'].astype('category')

        # Add patient info for nan samples
        # patient column resembls SAMPLE column but with numbers from 1 to 40
        m_patients = np.isnan(unpp_filtered_adata.obs['patient'])
        max_patient_old = unpp_filtered_adata.obs['patient'].max() + 1
        new_patients_anonym = np.unique(unpp_filtered_adata.obs.loc[m_patients, 'SAMPLE'])
        unpp_filtered_adata.obs.loc[m_patients, 'patient'] = unpp_filtered_adata.obs.loc[m_patients, 'SAMPLE']
        array_new_patients = np.arange(max_patient_old, max_patient_old + len(new_patients_anonym), 1, dtype=int)
        # 22 new patients: 78902 spots with 20613 genes
        # rename to 1-40
        unpp_filtered_adata.obs = unpp_filtered_adata.obs.replace({
            'patient': dict(zip(new_patients_anonym, array_new_patients))}, regex=True)
        unpp_filtered_adata.obs['patient'] = unpp_filtered_adata.obs['patient'].astype(int).astype('category')

        # Add batch == duplicates info same patient but different capture area
        # old samples replicates: same patient on same object slide
        # combination from lesional / non-lesional with patient ID
        # new samples replicates: same patient on multiple capture areas or object slides
        # combination lesional / non-lesional with patient ID
        unpp_filtered_adata.obs['replicates'] = 'Unknown_replicate_pairs'
        unpp_filtered_adata.obs['replicates'] = unpp_filtered_adata.obs[
            ['patient', 'biopsy_type']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)

        # to get the actual sample id == specimen: combine capture area with patient ID
        unpp_filtered_adata.obs['specimen'] = unpp_filtered_adata.obs['sample'].astype(str)
        unpp_filtered_adata.obs['specimen'] = unpp_filtered_adata.obs[
            ['capture_area', 'patient']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
        unpp_filtered_adata.obs['specimen'] = unpp_filtered_adata.obs['specimen'].astype('category')

        # save adata
        sc.write(os.path.join(adata_savepath, unppadata_cleaned_filename), unpp_filtered_adata)
        # unpp_filtered_adata.obs.to_excel(os.path.join(adata_savepath, 'MetaData_P15509_P16357_P21093.xlsx'))
    else:
        unpp_filtered_adata = sc.read(os.path.join(adata_savepath, unppadata_cleaned_filename))

    print("-------- Finished: Read out values --------")
    pp_adata, filename_adata = main(configs=configs_file, adata=unpp_filtered_adata, save_folder=output_path)
