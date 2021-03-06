from python_scripts.utils import helper_tools as ht

from operator import itemgetter
import numpy as np
import pandas as pd
import os
import sys


def get_spots_per_condition(adata, observable, cell_label, save_folder, key, paper_figure):
    """Read out spots for DGE analysis which fulfill specific conditions and save them in a .csv file

    Parameters
    ----------
    adata : annData
    observable : str
        observable in adata object
    cell_label : str
        skin layers (spatial transcriptomics data) or cell type label (single cell data)
    save_folder : str
    key : str
        if sc or spatial dataset
    paper_figure : str

    Returns
    -------

    """
    if "SC" in key:
        df_batches = dict()
        df_patients = dict()
        batches = adata.obs['sample'].values.categories
        for ind_batch, batch in enumerate(batches):
            # create dataframes
            df_batches[batch] = pd.DataFrame(columns=["sample", "batch"])
            df_patients[batch] = pd.DataFrame(columns=["sample", "batch"])

            # fill the dataframes
            df_batches[batch]['sample'] = adata.obs['sample'][adata.obs['sample'] == batch].values
            df_batches[batch]['batch'] = ind_batch + 1

            df_patients[batch]['sample'] = adata.obs['sample'][adata.obs['sample'] == batch].values
            df_patients[batch]['batch'] = 1
    else:
        # get batches assigned to each patient (currently we have 4 samples per patient)
        df_patients = ht.map_sample_batch_list(
            adata=adata, num_samples_patient=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 4])
        df_batches = ht.map_sample_batch_list(
            adata=adata, num_samples_patient=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4])

        # 1. Get all manual annotations by extracting all keys with upper case characters
        annotations = [char for char in adata.obs_keys() if any(c.isupper() for c in char)]

        if "ANNOTATOR" in annotations:
            annotations.remove("ANNOTATOR")

        elements = ["LESIONAL", "LESONAL", "NON LESIONAL"]
        index = []
        for el in elements:
            if el in annotations:
                index.append(annotations.index(el))
        try:
            disease_index = np.amin(index)
        except ValueError:
            print("REMINDER: Insert 'LESONAL' and 'NON LESIONAL' in your excel file")
            sys.exit()
        disease_labels = np.array(annotations[:disease_index])
        lesion_labels = np.array(itemgetter(*index)(annotations))

    # dataframe for DGE analysis
    df_condition = pd.DataFrame({"geneNames": adata.var.index})
    # meta Data
    df_metadata = pd.DataFrame(columns=["sample_id", "barcode", "condition", "label", "patient", "disease",
                                        "biopsy_type", "No_spots", "sizefactor", 'batch'])

    # spots having all required conditions fulfilled
    sample_names = np.unique(adata.obs['sample'])

    num_clusters = np.unique(adata.obs[observable])
    for n_clusters in num_clusters:
        # get count from spots and add counts to dataframe
        for ind, sample in enumerate(sample_names):
            # get a sub-dataframe for only the current sample
            ad = adata[adata.obs['sample'] == sample]
            if 'counts' in adata.layers.keys():
                sample_count_matrix = ad.layers['counts']
            else:
                # Read in normalised counts
                sample_count_matrix = ad.X

            # get index of spots in a cluster
            indexlist_obs_cluster = np.where(ad.obs[observable] == n_clusters)[0]

            # Important information for DGE Analysis
            no_spots = len(indexlist_obs_cluster)

            # save to metaData
            # "sample_id", "condition", "label" =(conditon_1 , ..), "batch", "disease", "No_spots"
            project = sample.split("_")[0]
            """ ATTENTION: This works only if we have one tissue biopsies per slide !! """
            if "SC" in key:
                disease_columns = np.array(['PSO'])
                type_lesional = np.array(['LESIONAL'])
                index_batch_sample = df_batches[project].index[df_batches[project]['sample'] == sample][0]
                index_patient_sample = df_patients[project].index[df_batches[project]['sample'] == sample][0]
            else:
                disease_columns = ad.obs[disease_labels][(ad.obs['sample'] == sample)] == 1
                disease_columns = np.asarray(disease_columns.idxmax(1).values)
                type_lesional = ad.obs[lesion_labels][(ad.obs['sample'] == sample)] == 1
                type_lesional = np.asarray(type_lesional.idxmax(1).values)
                index_batch_sample = df_batches[project].index[df_batches[project]['sample'] == sample][0]
                index_patient_sample = df_patients[project].index[df_patients[project]['sample'] == sample][0]

            # read out bulk vector of spots and label them in metaData with condition name
            # Single cell/spot Approach
            sample_barcodes = ad.obs.index
            matrix = sample_count_matrix[indexlist_obs_cluster]
            for index, spot in enumerate(matrix):
                df_condition[".".join(
                    [sample, sample_barcodes[index], observable, n_clusters]).replace(" ", "_")] = spot

                # save to metaData
                df_temp = pd.DataFrame({"sample_id": [sample],  # CaptureArea on slide
                                        "barcode": [sample_barcodes[index]],
                                        "condition": [n_clusters],
                                        "label": [ad.obs[cell_label].values[indexlist_obs_cluster[index]]],
                                        "patient": [int(df_patients[project]['batch'][index_patient_sample])],
                                        "disease": np.unique(disease_columns),
                                        "biopsy_type": np.unique(type_lesional),
                                        "No_spots": [no_spots],
                                        "sizefactor": [ad.obs['size_factors'].values[indexlist_obs_cluster[index]]],
                                        "batch": [int(df_batches[project]['batch'][index_batch_sample])]})  # Slide
                df_metadata = df_metadata.append(df_temp, ignore_index=True)

    df_metadata.to_csv(os.path.join(save_folder, "metaData_{}_{}.csv".format(observable, paper_figure)))
    df_condition.to_csv(os.path.join(save_folder, "{}_{}.csv".format(observable, paper_figure)))


def get_spots_per_condition_multiple(adata, observable, cell_label, save_folder, paper_figure):
    """Read out spots for DGE analysis which fulfill specific conditions and save them in a .csv file
        This function works also with multiple samples per capture area

    Parameters
    ----------
    adata : annData
    observable : str
        observable in adata object
    cell_label : str
        skin layers (spatial transcriptomics data) or cell type label (single cell data)
    save_folder : str
    paper_figure : str

    Returns
    -------

    """
    # dataframe for DGE analysis
    df_condition = pd.DataFrame({"geneNames": adata.var.index})
    # meta Data
    df_metadata = pd.DataFrame(
        columns=["barcode", "sample", "project", "specimen", "condition", "cluster_label", "patient", "disease",
                 "biopsy_type", "n_spots", "sizefactor", 'batch', 'capture_area'])

    # spots having all required conditions fulfilled
    specimen_names = np.unique(adata.obs['specimen'])

    cluster_labels = np.unique(adata.obs[observable])
    for cluster_label in cluster_labels:
        # get count from spots and add counts to dataframe
        for ind, specimen in enumerate(specimen_names):
            # get a sub-dataframe for only the current specimen
            ad = adata[adata.obs['specimen'] == specimen].copy()
            if 'counts' in adata.layers.keys():
                sample_count_matrix = ad.layers['counts']
            else:
                # Read in normalised counts
                sample_count_matrix = ad.X

            # get index of spots in a cluster
            indexlist_obs_cluster = np.where(ad.obs[observable] == cluster_label)[0]

            # Important information for DGE Analysis
            n_spots = len(indexlist_obs_cluster)

            # read out bulk vector of spots and label them in metaData with condition name
            # Single cell/spot Approach
            sample_barcodes = ad.obs.index
            matrix = sample_count_matrix[indexlist_obs_cluster]
            for index, spot in enumerate(matrix):
                # Must be named like that as index in obs is not unique, only unique on a single capture area
                df_condition[".".join(
                    [specimen, sample_barcodes[index], observable, cluster_label]).replace(" ", "_")] = spot

                # save to metaData
                df_temp = pd.DataFrame({
                    "barcode": [sample_barcodes[index]],
                    "sample": [ad.obs['sample'].values[indexlist_obs_cluster[index]]],
                    "project": [ad.obs['project'].values[indexlist_obs_cluster[index]]],
                    "specimen": [specimen],
                    "condition": [cluster_label],
                    "cluster_label": [ad.obs[cell_label].values[indexlist_obs_cluster[index]]],
                    "patient": [int(ad.obs['patient'].values[indexlist_obs_cluster[index]])],
                    "disease": [ad.obs['DISEASE'].values[indexlist_obs_cluster[index]]],
                    "biopsy_type": [ad.obs['biopsy_type'].values[indexlist_obs_cluster[index]]],
                    "n_spots": [n_spots],
                    "sizefactor": [ad.obs['size_factors'].values[indexlist_obs_cluster[index]]],
                    "batch": [ad.obs['object_slide'].values[indexlist_obs_cluster[index]]],  # Slide
                    "capture_area": [ad.obs['capture_area'].values[indexlist_obs_cluster[index]]]})  # CaptureArea on slide
                df_metadata = df_metadata.append(df_temp, ignore_index=True)

    df_metadata.to_csv(os.path.join(save_folder, "metaData_{}_{}.csv".format(observable, paper_figure)))
    df_condition.to_csv(os.path.join(save_folder, "{}_{}.csv".format(observable, paper_figure)))
