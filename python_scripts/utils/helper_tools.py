import os
import json
import numpy as np
import pandas as pd
from operator import itemgetter
import random
import copy
import math
import configparser
import glob


def load_config(config_path):
    # Load the configuration file
    configfile_list = glob.glob(os.path.join(config_path, "*.ini"))
    if len(configfile_list) > 1:
        raise ValueError("More than one config file found: {}".format(configfile_list))

    with open(configfile_list[0]) as f:
        sample_config = f.read()
    config_file = configparser.RawConfigParser(allow_no_value=True)
    config_file.read_string(sample_config)

    return config_file


def get_files(path, ending):
    list_csv_files = []

    for file in os.listdir(path):
        if file.endswith(ending):
            list_csv_files.append(os.path.join(path, file))

    return list_csv_files


def get_list(dictionary):
    """
    From dictionary get keys as list

    :param dictionary: dict
    :return: list
    """

    key_list = []
    for key in dictionary.keys():
        key_list.append(key)

    return key_list


def get_indexes(dfobj, value):
    """
    Get index positions of value in dataframe i.e. dfObj.

    :param dfobj:
    :param value:
    :return:
    """

    listofpos = list()
    # Get bool dataframe with True at positions where the given value exists
    result = dfobj.isin([value])
    # Get list of columns that contains the value
    seriesobj = result.any()
    columnnames = list(seriesobj[seriesobj is True].index)
    # Iterate over list of columns and fetch the rows indexes where value exists
    for col in columnnames:
        rows = list(result[col][result[col] is True].index)
        for row in rows:
            listofpos.append((row, col))

    # Return a list of tuples indicating the positions of value in the dataframe
    return listofpos


def load_sample_config_file(filename, file_type):
    """
    read in configs like paths, files names

    :param filename:
    :param file_type:
    :return: read configs
    """

    if file_type is "csv":
        configs = pd.read_csv(filename, sep=";")

        # No need to save configs in different variables; it is just for visualization and easier to read out later
        file_path = configs['file_path'].values  # 0
        library_path = configs['output_type_spaceranger'].values  # 1
        sample_strings = configs['sample_strings'].values  # 2
        sample_id_strings = configs['sample_id_strings'].astype(str).values  # 3

        # output matrix
        output_type_matrix = configs['output_type_matrix'].values  # 4
        raw_feature_bc_matrix_folder = configs['raw_feature_bc_matrix_folder'].values  # 5
        filtered_feature_bc_matrix_folder = configs['filtered_feature_bc_matrix_folder'].values  # 6

        barcode_file_end = configs['barcode_file_end'].values  # 7
        features_file_end = configs['features_file_end'].values  # 8
        matrix_file_end = configs['matrix_file_end'].values  # 9
        raw_hdf5_file_end = configs['raw_hdf5_file_end'].values  # 10
        filtered_hdf5_file_end = configs['filtered_hdf5_file_end'].values  # 11

        # spatial information
        output_type_spatial = configs['output_type_spatial'].values  # 12

        tissue_pos_list_file_end = configs['tissue_pos_list_file_end'].values  # 13
        scale_factors_file_end = configs['scale_factors_file_end'].values  # 14
        aligned_fiducials_file_end = configs['aligned_fiducials_file_end'].values  # 15
        detected_tissue_file_end = configs['detected_tissue_file_end'].values  # 16
        tissue_hires_file_end = configs['tissue_hires_file_end'].values  # 17
        tissue_low_res_file_end = configs['tissue_low_res_file_end'].values  # 18

        # annotation files
        annotation_path = configs["annotation_path"].values  # 19
        excel_sheet = configs['excel_sheet'].values  # 20

        # lesional_folder = configs['lesional_folder'].values  # 21
        # non_lesional_folder = configs['non-lesional_folder'].values  # 22

        # velocyto (loom files)
        velocyto_path = configs['velocity_path'].values  # 23 -> 21
        loom_file = configs['loom_file_end'].values  # 24 -> 22

        # Capture area
        capture_area = configs['capture_area'].values  # 25 -> 23
        patient = configs['patient'].values  # 26 -> 24

    else:
        with open(filename) as json_file:
            configs = json.load(json_file)

            # No need to save configs in different variables; it is just for visualization and easier to read out later
            file_path = configs['file_path']  # 0
            library_path = configs['output_type_spaceranger']  # 1
            sample_strings = configs['sample_strings']  # 2
            sample_id_strings = configs['sample_id_strings']  # 3

            # output matrix
            output_type_matrix = configs['output_type_matrix']  # 4
            raw_feature_bc_matrix_folder = configs['raw_feature_bc_matrix_folder']  # 5
            filtered_feature_bc_matrix_folder = configs['filtered_feature_bc_matrix_folder']  # 6

            barcode_file_end = configs['barcode_file_end']  # 7
            features_file_end = configs['features_file_end']  # 8
            matrix_file_end = configs['matrix_file_end']  # 9
            raw_hdf5_file_end = configs['raw_hdf5_file_end']  # 10
            filtered_hdf5_file_end = configs['filtered_hdf5_file_end']  # 11

            # spatial information
            output_type_spatial = configs['output_type_spatial']  # 12

            tissue_pos_list_file_end = configs['tissue_pos_list_file_end']  # 13
            scale_factors_file_end = configs['scale_factors_file_end']  # 14
            aligned_fiducials_file_end = configs['aligned_fiducials_file_end']  # 15
            detected_tissue_file_end = configs['detected_tissue_file_end']  # 16
            tissue_hires_file_end = configs['tissue_hires_file_end']  # 17
            tissue_low_res_file_end = configs['tissue_low_res_file_end']  # 18

            # annotation files
            annotation_path = configs["annotation_path"]  # 19
            excel_sheet = configs['excel_sheet']  # 20

            # lesional_folder = configs['lesional_folder']  # 21
            # non_lesional_folder = configs['non-lesional_folder']  # 22

            # velocyto (loom files)
            velocyto_path = configs['velocity_path']  # 23
            loom_file = configs['loom_file_end']  # 24

            # Capture area
            capture_area = configs['capture_area'].values  # 25
            patient = configs['patient'].values  # 26

    return \
        file_path, library_path, sample_strings, sample_id_strings, output_type_matrix, \
        raw_feature_bc_matrix_folder, filtered_feature_bc_matrix_folder, barcode_file_end, features_file_end, \
        matrix_file_end, raw_hdf5_file_end, filtered_hdf5_file_end, output_type_spatial, tissue_pos_list_file_end, \
        scale_factors_file_end, aligned_fiducials_file_end, detected_tissue_file_end, tissue_hires_file_end, \
        tissue_low_res_file_end, annotation_path, excel_sheet, velocyto_path, loom_file, capture_area, patient


def save_obsm_with_library_id(adatas):
    # save .obsm (image coordiantes) as list containing [['X' x 'Y' x 'sample_id'], ...]
    keys_obsm = []
    if len(adatas) > 2:
        list_adata_obsm = adatas[0].obsm['spatial']
        for obj in adatas[1:]:
            unique_samples = np.unique(obj.obs['sample'])
            for sample_name in unique_samples:
                # key_list = obj.obsm_keys()[0]
                size_obj = obj.obsm['spatial'].shape[0]
                keys_obsm.append(sample_name)
                lib_id = np.array(size_obj * [sample_name])
                # convert coordinates to string and save coordinates with library_id
                str_coords = [[str(var[0]), str(var[1]), lib_id[ind]] for ind, var in enumerate(obj.obsm['spatial'])]
                list_adata_obsm.extend(str_coords)
    else:
        list_adata_obsm = []
        for obj in adatas:
            unique_sample = np.unique(obj.obs['sample'])[0]
            # key_list = obj.obsm_keys()[0]
            size_obj = obj.obsm['spatial'].shape[0]
            keys_obsm.append(unique_sample)
            lib_id = np.array(size_obj * [unique_sample])
            # convert coordinates to string and save coordinates with library_id
            str_coords = [[str(var[0]), str(var[1]), lib_id[ind]] for ind, var in enumerate(obj.obsm['spatial'])]
            list_adata_obsm.extend(str_coords)


def concatenate_adatas(adatas, batch_key='library_id'):
    """
    Concatenate multiple adatas by library_id
    CARE: adata.uns['spatial'][library_id] must be unique for each adata! otherwise same properties will be overwritten
    ==> save .uns and .obsm and before and then load it into concatenated annData object

    :param adatas: [list] of adata objects
    :param batch_key: observable name to add to distinguish between batches
    :return: return single adata object containing all information including images
    """
    adata_spatial = adatas[0].copy()

    # save .obsm (image coordiantes) as list containing [[X x Y], ...]
    # to concatenate loom files also save .layers
    list_adata_obsm = []
    for obj in adatas:
        list_adata_obsm.extend(obj.obsm['spatial'])

    # save .uns (images, image properties)
    list_adata_uns = [ann_data.uns['spatial'] for ann_data in adatas]
    adata_spatial = adata_spatial.concatenate(adatas[1:], batch_key=batch_key,
                                              batch_categories=[k + '_' + str(counter + 1)
                                                                for counter, d in enumerate(list_adata_uns)
                                                                for k, v in d.items()],
                                              uns_merge='unique')

    # obsm: now get coordinates by checking third entry of adata_spatial.obsm['spatial'] and then convert strings to int
    # adata_spatial.obsm['spatial'] = np.asarray(list_adata_obsm)

    return adata_spatial


def store_categories_as_clusters(adata, annotations):
    """
    Determine one label per spot, by choosing always the category with lowest No. of spots assigned to it.
    Save it in a vector and add it as new observable to annData object.

    :param adata: [annData]
    :param annotations: [list]
    :return: [annData]
    """
    # split annotations by disease and tissue / cell types
    elements = ["LESONAL", "NON LESIONAL"]
    try:
        index = []
        for el in elements:
            index.append(annotations.index(el))
        target_index = np.amax(index) + 1
    except ValueError:
        target_index = None
        print("REMINDER: Insert 'LESONAL' and 'NON LESIONAL' in your excel file")

    tissue_cell_labels = annotations[target_index:]

    # remove cell cycle annotations from tissue_cell_labels list
    for cc in ["G1_score", "G2M_score", "S_score", "M_score"]:
        try:
            tissue_cell_labels.remove(cc)
        except ValueError:
            continue

    # save annotations in one vector
    # if one spot got several labels choose the one which has the lowest number of total spots assigned
    if len(tissue_cell_labels) > 0:
        if "ANNOTATOR" in tissue_cell_labels:
            tissue_cell_labels.remove("ANNOTATOR")
        spots_per_label = dict()
        for label in tissue_cell_labels:
            spots_per_label[label] = len(adata.obs[label][adata.obs[label] == 1])

        value_missing_assignments = len(tissue_cell_labels) + 1
        label_missing_assignments = "Unknown"
        label_per_spot = []
        value_per_spot = []
        # loop through each spot (super slow method ...)
        for c in range(len(adata.obs.index.values)):
            # get labels from that spot and label it with the category having the lowest number of spots assigned to it
            spot_adata = adata[c]
            # get non zero indices and label name
            indices = np.nonzero(spot_adata.obs[tissue_cell_labels].to_numpy())[1]

            if len(indices) > 1:
                labels = list(itemgetter(*indices)(tissue_cell_labels))

                # TODO find a better way to choose between multiple labels the one which shall be used for the spot
                sequence_selection = np.array([i for i in range(len(labels))])
                # define weights using the number of spots assigned to these labels
                spots_counts = [spots_per_label.get(cat) for cat in labels]
                index_spot = random.choices(sequence_selection, weights=1/np.asarray(spots_counts), k=1)[0]

                # save cluster label with int values from 0 to number of tissue/cell type annotations
                value_per_spot.append(indices[index_spot])
                label_per_spot.append(labels[index_spot])
            elif len(indices) == 1:
                labels = tissue_cell_labels[indices[0]]
                # save cluster label with int values from 0 to number of tissue/cell type annotations
                value_per_spot.append(indices[0])
                label_per_spot.append(labels)
            else:
                # spots without label ...
                label_per_spot.append(label_missing_assignments)
                value_per_spot.append(value_missing_assignments)

        # missing spots barcodes index
        missing_spots_index = np.where(np.asarray(label_per_spot) == label_missing_assignments)[0]
        if len(missing_spots_index) > 0:
            missing_spots_samples = np.unique(adata.obs['sample'][
                                                  np.where(np.asarray(label_per_spot) == label_missing_assignments)[0]])
            print("Missing spots from sample(s): ", missing_spots_samples)

        # save labels per spot as observable in adata object
        adata.obs['groups'] = label_per_spot
        adata.obs['cat_clusters'] = value_per_spot
        # store lables as categories
        adata.obs['groups'] = adata.obs['groups'].astype('category')
        adata.obs['cat_clusters'] = adata.obs['cat_clusters'].astype('category')

    else:
        print("No manual annotations found")

    return adata


def map_sample_batch_list(adata, num_samples_patient):
    """

    :param adata:
    :param num_samples_patient: [list] [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2]
    :return:
    """

    def order_of_magnitude(number):
        return math.floor(math.log(number, 10))

    metadata_batch = pd.DataFrame(columns=["sample", "batch"])
    metadata = dict()
    list_project_samples = dict()
    counter_batch = 1
    for ind, project in enumerate(np.unique(adata.obs['project'])):
        list_project_samples[project] = np.unique(adata.obs['sample'][adata.obs['project'] == project])
        sample_names = [int(i.split('_')[1]) for i in list_project_samples[project]]

        temp_num_oom = 10 ** order_of_magnitude(np.amax(sample_names))
        max_samplenumber = np.amax(sample_names) - temp_num_oom

        rest = max_samplenumber % num_samples_patient[ind]
        if rest != 0:
            max_samplenumber = max_samplenumber + rest

        # determine no batches
        num_batches = 0
        list_batches = num_samples_patient.copy()
        number_samples = copy.copy(max_samplenumber)
        while True:
            list_batches.remove(num_samples_patient[num_batches])
            if ((number_samples - num_samples_patient[num_batches]) == 0) | \
                    (num_batches == (len(num_samples_patient) - 1)):
                num_batches += 1
                break
            number_samples = number_samples - num_samples_patient[num_batches]
            num_batches += 1

        batches = []
        for i_c in range(num_batches):
            batches.extend([i_c + 1 + ind] * num_samples_patient[i_c])

        array_sampleids = np.arange(1 + temp_num_oom, temp_num_oom + max_samplenumber + 1)

        for i_c in range(len(array_sampleids)):
            metadata_batch.loc[i_c] = ["_".join([project, str(array_sampleids[i_c])]), batches[i_c]]

        metadata[project] = pd.DataFrame(columns=["sample", "batch"])
        for i_c in range(len(sample_names)):
            sample_batch = metadata_batch[(metadata_batch ==
                                           list_project_samples[project][i_c]).any(1)].stack().unique()
            metadata[project].loc[i_c] = [sample_batch[0], sample_batch[1]]
        num_samples_patient = list_batches

    return metadata


def get_tissue_annot(adata):
    # Use annotations from pathologist instead of clusters
    obs_keys = list(adata.obs_keys())
    # get all manual annotations by extracting all keys with upper case characters
    annotations = [char for char in obs_keys if any(c.isupper() for c in char)]

    if "ANNOTATOR" in annotations:
        annotations.remove("ANNOTATOR")

    # split annotations by disease and tissue / cell types
    elements = ["LESONAL", "NON LESIONAL", "LESIONAL"]
    intersection_ele = np.intersect1d(adata.obs.columns, elements)
    try:
        index = []
        for el in intersection_ele:
            index.append(annotations.index(el))
        target_index = np.amax(index) + 1
        disease_index = np.amin(index)
    except ValueError:
        target_index = None
        print("REMINDER: Insert 'LESONAL' and 'NON LESIONAL' in your excel file")

    tissue_cell_labels = annotations[target_index:]
    disease_labels = np.array(annotations[:disease_index])
    lesion_labels = np.array(itemgetter(*index)(annotations))

    # remove cell cycle annotations from tissue_cell_labels list
    for cc in ["G1_score", "G2M_score", "S_score", "M_score"]:
        try:
            tissue_cell_labels.remove(cc)
        except ValueError:
            continue

    return tissue_cell_labels, disease_labels, lesion_labels


def assign_donor_spot(adata, no_samples_per_patient=4):
    """
    Only usable if Number of samples per patient is equal across whole dataset

    :param adata:
    :param no_samples_per_patient:
    :return:
    """
    # get number of samples in data set and divide it by 4 (because each patient donored 4 samples)
    samples = np.unique(adata.obs['library_id'])
    no_donors = int(len(samples) / no_samples_per_patient)
    list_donors = np.arange(1, no_donors + 1).astype(str).tolist()
    adata.obs['donor'] = adata.obs['library_id'].cat.add_categories(list_donors)
    counter = 0
    for c_donors in range(0, len(samples), no_samples_per_patient):
        counter += 1
        donor_samples = samples[c_donors: no_samples_per_patient + c_donors]
        adata.obs['donor'][np.in1d(adata.obs['library_id'], donor_samples)] = str(counter)
    adata.obs['donor'] = adata.obs['donor'].cat.remove_unused_categories()
