#! /usr/bin/python
"""Create config file with used parameters for analysis
    File name: init_variables.py
    Author: Christina Hillig
    Date created: April/29/2021
    Date last modified: May/02/2021
    Python Version: 3.7
"""

from datetime import date
import os
import sys
import numpy as np
try:
    import configparser  # for python3
except ImportError:
    import ConfigParser as configparser  # for python2


def create_config(path):
    """Create config file
    Always save and load config with adata object

    Parameters
    ----------
    path : str

    Returns
    -------

    """
    configfile_name = os.path.join(path, "P15509_P16357_config.ini")

    # Check if there is already a configuration file
    if not os.path.isfile(configfile_name):
        # Create the configuration file as it doesn't exist yet
        cfgfile = open(configfile_name, "w")

        # Add content to the file
        config = configparser.ConfigParser()
        config.add_section("project")
        config.set("project", "project_name", "Immune Publication")
        config.set("project", "author", "Christina Hillig")
        config.set("project", "date", str(date.today()))

        config.add_section("data")
        config.set("data", "data_type", "Spatial Transcriptomics")
        # Path to Database
        database_path = \
            'Volumes' + os.sep + 'Samsung_T5' + os.sep + 'IGSSE-TUM__datasets' + os.sep + \
            '2020__Kilian_Eyerich__Spatial_transcriptomics_skin_single-cell_RNA-seq' + os.sep + 'Database'
        config.set("data", "input_path", database_path)
        config.set("data", "output_path", path)
        config.set("data", "sequencing_samples", '["P15509", "P16357"]')
        # 1. ids for project number P16357
        ids = np.arange(1001, 1061, 1)
        # 2. all ids from project P15509 and P16357
        sample_id = ['1001', '1002', '1003', '1004'] + list(ids.astype(str))
        config.set("data", "sample_ids", str(sample_id))

        config.add_section("preprocessing")
        config.set("preprocessing", "read_raw_matrix", 'False')
        # combined analysis on all samples if from different data sets and samples separated into batches
        config.set("preprocessing", "sample_concat", 'True')
        config.set("preprocessing", "apply_mt_threshold", 'True')
        # Set to True if you work with single cell data
        config.set('preprocessing', 'filter_doublets', 'False')
        # which kind of normalisation step shall be applied: scran (R) or scanpy
        config.set("preprocessing", "normalisation_function", 'scran')
        config.set("preprocessing", "exclude_highly_expressed", 'False')
        config.set("preprocessing", "apply_scaling", 'False')
        config.set("preprocessing", "get_cc_effect", 'False')
        config.set("preprocessing", "apply_remove_cc_effect", 'False')
        config.set("preprocessing", "apply_batch_correction", "True")
        config.set("preprocessing", "batch_correction_tool", "scanorama")

        config.add_section("analysis")

        config.add_section("input_files")
        config.set("input_files", "cell_cycle",
                   os.path.join("..", "..", 'input', 'cell_cycle', 'Macosko_cell_cycle_genes_2015.csv'))
        config.set("input_files", "rawdata_paths",
                   os.path.join("..", "..", 'input', 'rawdata_paths', 'P15509_P16357_wo_4_7.csv'))

        config.add_section("visualisation_options")
        config.set("visualisation_options", "image_res", "hires")

        # config.add_section("Set_parameters")

        # Write config file and save with adata
        config.write(cfgfile)
        cfgfile.close()


def create_config_review(path):
    """Create config file
    Always save and load config with adata object

    Parameters
    ----------
    path : str

    Returns
    -------

    """
    configfile_name = os.path.join(path, "P15509_P16357_P21093_config.ini")

    # Check if there is already a configuration file
    if not os.path.isfile(configfile_name):
        # Create the configuration file as it doesn't exist yet
        cfgfile = open(configfile_name, "w")

        # Add content to the file
        config = configparser.ConfigParser()
        config.add_section("project")
        config.set("project", "project_name", "Immune Publication")
        config.set("project", "author", "Christina Hillig")
        config.set("project", "date", str(date.today()))

        config.add_section("data")
        config.set("data", "data_type", "Spatial Transcriptomics")
        # Path to Database
        database_path = \
            'Volumes' + os.sep + 'Samsung_T5' + os.sep + 'IGSSE-TUM__datasets' + os.sep + \
            '2020__Kilian_Eyerich__Spatial_transcriptomics_skin_single-cell_RNA-seq' + os.sep + 'Database'
        config.set("data", "input_path", database_path)
        config.set("data", "output_path", path)
        config.set("data", "sequencing_samples", '["P15509", "P16357", "P21093"]')
        # 1. ids for project number P16357
        ids = np.arange(1001, 1061, 1)
        # 2. all ids from project P15509 and P16357
        sample_id = ['1001', '1002', '1003', '1004'] + list(ids.astype(str))
        # 3. ids from project P15509 and P16357 and P21093
        ids_p2 = np.arange(8958, 8977, 1)
        ids_p2 = [str(s) + '21L00' for s in ids_p2]
        sample_id = sample_id + ids_p2
        config.set("data", "sample_ids", str(sample_id))

        config.add_section("preprocessing")
        config.set("preprocessing", "read_raw_matrix", 'False')
        # combined analysis on all samples if from different data sets and samples separated into batches
        config.set("preprocessing", "sample_concat", 'True')
        config.set("preprocessing", "apply_mt_threshold", 'True')
        # Set to True if you work with single cell data
        config.set('preprocessing', 'filter_doublets', 'False')
        # which kind of normalisation step shall be applied: scran (R) or scanpy
        config.set("preprocessing", "normalisation_function", 'scran')
        config.set("preprocessing", "normalisation_groups", 'tissue_layer')
        config.set("preprocessing", "exclude_highly_expressed", 'False')
        config.set("preprocessing", "apply_scaling", 'False')
        config.set("preprocessing", "get_cc_effect", 'False')
        config.set("preprocessing", "apply_remove_cc_effect", 'False')
        config.set("preprocessing", "apply_batch_correction", "True")
        config.set("preprocessing", "batch_correction_tool", "scanorama")

        config.add_section("analysis")

        config.add_section("input_files")
        config.set("input_files", "cell_cycle",
                   os.path.join("..", "..", 'input', 'cell_cycle', 'Macosko_cell_cycle_genes_2015.csv'))
        config.set("input_files", "rawdata_paths",
                   os.path.join("..", "..", 'input', 'rawdata_paths', 'P15509_P16357_P21093_wo_4_7.csv'))

        config.add_section("visualisation_options")
        config.set("visualisation_options", "image_res", "hires")

        # config.add_section("Set_parameters")

        # Write config file and save with adata
        config.write(cfgfile)
        cfgfile.close()


def init_vars():
    # initialization variables
    # create saving folder in current project path
    # today = date.today()
    today = '2022-04-08'
    savepath = os.path.join("..", "..", "adata_storage", str(today))
    os.makedirs(savepath, exist_ok=True)

    # Create config file
    create_config_review(path=savepath)

    return savepath


if __name__ == '__main__':
    init_vars()


