from python_scripts.utils import helper_tools as ht

# packages
import h5py
import json
import numpy as np
import pandas as pd
import os
import csv
import collections
import tables
import scipy.io
import scipy.sparse as sp_sparse
from tqdm import tqdm

# Single cell packages
import scanpy as sc
from scanpy import logging as logg
import scvelo as scv

# Plot packages
from matplotlib.image import imread


def _read_excel(adata, path_excel_file):
    """
    Load Kilians and Alex annotation files
    keys :  'SPOT', 'ANNOTATOR'
            'PSO', 'AE', 'LICHEN', 'PRP', 'LESONAL', 'NON LESIONAL',
            'upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS', 'DERMIS',
            'DERdepth1', 'DERdepth2', 'DERdepth3', 'DERdepth4', 'DERdepth5', 'DERdepth6', 'DERdepth7',
            'INTERFACE', 'VESSEL', 'HAIRFOLLICLE', 'GLAND', 'NEUTROPHIL', 'MONONUCLEAR', 'MUSCLE', 'KERATINOCYTE',
            'FIBROBLAST', 'ENDOTHELIAL', 'CONNECTIVE TISSUE', 'MELANOPHAGE'

    :param adata: [annData]
    :param path_excel_file: [string]
    :return: [annData]
    """

    annotation_file = pd.read_excel(path_excel_file)

    annotation_file.rename(columns={'SPOT': 'barcode'}, inplace=True)
    annotation_file.set_index('barcode', inplace=True)

    adata.obs = adata.obs.join(annotation_file, how="right")

    return adata


def _read_visium(adata, path, library_id=None, load_images=True):
    """
    Read visium data like done in scanpy's sc.read_visium() but with dense annData object

    :param adata: [annData]
    :param path: [string]
    :param library_id: [string]
    :param load_images: [bool]
    :return: [annData]
    """

    adata.uns["spatial"] = dict()

    from h5py import File

    with File(path[0], mode="r") as f:
        attrs = dict(f.attrs)

    if library_id is None:
        library_id = str(attrs.pop("library_ids")[0], "utf-8")

    adata.uns["spatial"][library_id] = dict()

    if load_images:
        files = dict(
            tissue_positions_file=os.path.join(path[1], 'tissue_positions_list.csv'),
            scalefactors_json_file=os.path.join(path[1], 'scalefactors_json.json'),
            hires_image=os.path.join(path[1], 'tissue_hires_image.png'),
            lowres_image=os.path.join(path[1], 'tissue_lowres_image.png'),
        )

        # check if files exists, continue if images are missing
        for f in files.values():
            if not os.path.exists(f):
                if any(x in str(f) for x in ["hires_image", "lowres_image"]):
                    logg.warning(
                        "You seem to be missing an image file.\n"
                        "Could not find '{}'.".format(f))
                else:
                    raise OSError("Could not find '{}'".format(f))

        adata.uns["spatial"][library_id]['images'] = dict()
        for res in ['hires', 'lowres']:
            try:
                adata.uns["spatial"][library_id]['images'][res] = imread(str(files['{}_image'.format(res)]))
            except Exception:
                raise OSError("Could not find '{}_image'".format(res))

        # read json scalefactors
        json_file_path = files['scalefactors_json_file']
        with open(json_file_path) as json_file:
            adata.uns["spatial"][library_id]['scalefactors'] = json.loads(json_file.read())

        adata.uns["spatial"][library_id]["metadata"] = {
            k: (str(attrs[k], "utf-8") if isinstance(attrs[k], bytes) else attrs[k])
            for k in ("chemistry_description", "software_version")
            if k in attrs
        }

        # read coordinates
        positions = pd.read_csv(files['tissue_positions_file'], header=None)
        positions.columns = [
            'barcode',
            'in_tissue',
            'array_row',
            'array_col',
            'pxl_col_in_fullres',
            'pxl_row_in_fullres',
        ]
        positions.index = positions['barcode']

        adata.obs = adata.obs.join(positions, how="left")

        adata.obsm['spatial'] = adata.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()

        adata.obs.drop(columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'], inplace=True, )

    return adata


def _get_matrix_from_h5(filename):
    """
    10x Genomics version
    https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/advanced/h5_matrices

    :param filename: [string]
    :return:
    """
    if filename is not None:

        CountMatrix = collections.namedtuple('CountMatrix', ['feature_ref', 'barcodes', 'matrix'])
        with tables.open_file(filename, 'r') as f:
            mat_group = f.get_node(f.root, 'matrix')
            barcodes = f.get_node(mat_group, 'barcodes').read()
            data = getattr(mat_group, 'data').read()
            indices = getattr(mat_group, 'indices').read()
            indptr = getattr(mat_group, 'indptr').read()
            shape = getattr(mat_group, 'shape').read()
            matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)

            feature_ref = {}
            feature_group = f.get_node(mat_group, 'features')
            feature_ids = getattr(feature_group, 'id').read()
            feature_names = getattr(feature_group, 'name').read()
            feature_types = getattr(feature_group, 'feature_type').read()
            feature_ref['id'] = feature_ids
            feature_ref['name'] = feature_names
            feature_ref['feature_type'] = feature_types
            tag_keys = getattr(feature_group, '_all_tag_keys').read()
            feature_ref['tag_keys'] = tag_keys
            # for key in tag_keys:
            #     feature_ref[key] = getattr(feature_group, key).read()

            return CountMatrix(feature_ref, barcodes, matrix)
    else:
        return None


def _get_feature_bc_matrix(matrix_dir):
    """
    10x Genomics version:
    https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices

    :param matrix_dir: [string]
    :return:
    """

    mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))

    features_path = os.path.join(matrix_dir, "features.tsv")
    feature_ids = [row[0] for row in csv.reader(features_path, delimiter="\t")]
    gene_names = [row[1] for row in csv.reader(features_path, delimiter="\t")]
    feature_types = [row[2] for row in csv.reader(features_path, delimiter="\t")]
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv")
    barcodes = [row[0] for row in csv.reader(barcodes_path, delimiter="\t")]

    return mat, feature_ids, gene_names, feature_types, barcodes


def _read_hdf5_files(path_h5_file):
    """
    # TODO WIP
    My version

    :param path_h5_file: [string]
    :return: [Dataframe] pandas dataframe
    """
    # path_h5_file = "/Volumes/Samsung_T5/IGSSE-TUM__datasets/" \
    #                "2020__Kilian_Eyerich__Spatial_transcriptomics_skin_single-cell_RNA-seq/" \
    #                "delivery_project_03061/fromBCF/P15509_skin/sr_rev0/spaceranger/P15509_1001/" \
    #                "count_matrices"
    h5_filtered_feature_bc_matrix = os.path.join(path_h5_file, "filtered_feature_bc_matrix.h5")

    list_parent_keys_h5_filtered = []
    list_child_keys_h5_filteres = {}
    matrix_ele = {}
    with h5py.File(h5_filtered_feature_bc_matrix, 'r') as f:
        # List all groups
        for key in f.keys():
            print(key)
            list_parent_keys_h5_filtered.append(key)

        print("Keys: %s" % f.keys())

        # Get the data list(f[list_parent_keys_h5_filtered[0]]['barcodes'])
        for c, var in enumerate(list_parent_keys_h5_filtered):
            list_child_keys_h5_filteres[var] = []
            list_child_keys_h5_filteres[var].append(list(f[var].keys()))
            for i, var_child in enumerate(list_child_keys_h5_filteres[var][0]):
                matrix_ele[var_child] = []
                matrix_ele[var_child].append(list(f[var][var_child]))

    # create csv file from h5 file
    spatialtranscriptomics_h5_data = {
        'barcodes': matrix_ele[list_child_keys_h5_filteres[list_parent_keys_h5_filtered[0]][0][0]][0],
        'data': matrix_ele[list_child_keys_h5_filteres[list_parent_keys_h5_filtered[0]][0][1]][0],
        'features': matrix_ele[list_child_keys_h5_filteres[list_parent_keys_h5_filtered[0]][0][2]][0],
        'indices': matrix_ele[list_child_keys_h5_filteres[list_parent_keys_h5_filtered[0]][0][3]][0],
        'indptr': matrix_ele[list_child_keys_h5_filteres[list_parent_keys_h5_filtered[0]][0][4]][0],
        'shape': matrix_ele[list_child_keys_h5_filteres[list_parent_keys_h5_filtered[0]][0][5]][0]}

    spatialtranscriptomics_h5_df = pd.DataFrame(spatialtranscriptomics_h5_data, columns=['barcodes', 'data',
                                                                                         'features', 'indices',
                                                                                         'indptr', 'shape'])
    # spatialtranscriptomics_h5_df.to_csv(os.path.join(path_h5_file, "spatialtranscriptomics_h5.csv"), sep="\t")

    return spatialtranscriptomics_h5_df


def _read_feature_bc_matrix_tsv(filename):
    """
    My version

    :param filename: [string]
    :return: [Dataframe] pandas data frame, Compressed Sparse Column Format (CSC)
    """
    # define file specific parts
    mat = scipy.io.mmread(os.path.join(filename, "matrix.mtx"))
    features_path = os.path.join(filename, "features.tsv")
    barcodes_path = os.path.join(filename, "barcodes.tsv")

    # read in files with csv package
    # pd.read_csv(features_path, delimiter="\t", header=None)
    feature_ids = pd.read_csv(features_path, delimiter="\t", header=None)[0]
    gene_names = pd.read_csv(features_path, delimiter="\t", header=None)[1]
    feature_types = pd.read_csv(features_path, delimiter="\t", header=None)[2]
    barcodes = pd.read_csv(barcodes_path, delimiter="\t", header=None)[0]

    # create dataFrame to save features in one csv file
    spatialtranscriptomics_data = {'feature_ids': feature_ids,
                                   'gene_names': gene_names,
                                   'feature_types': feature_types,
                                   'barcodes': barcodes}

    spatialtranscriptomics_df = pd.DataFrame(spatialtranscriptomics_data, columns=['feature_ids', 'gene_names',
                                                                                   'feature_types', 'barcodes'])
    # save barcodes and features to csv file
    spatialtranscriptomics_df.to_csv(os.path.join(filename,
                                                  "spatialtranscriptomics_countmatrix.csv"), sep="\t")

    return spatialtranscriptomics_data, mat


def _scanpy_load_annotate_tsv_mtx_files(path_filtered_files, filenames, path_poslist_file, path_raw_files=None,
                                        read_raw_matrix=False):
    """
    source: https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/
    Case-study_Mouse-intestinal-epithelium_1906.ipynb

    :param path_raw_files: [string] contains path to raw mtx and tsv files created by 10x Genomics Spaceranger
    :param path_filtered_files: [string] contains path to filtered mtx and tsv files created by 10x Genomics Spaceranger
    :param filenames: [string] contains the names of the mtx and tsv files in the following order:
    [matrix.mtx, features.tsv, barcodes.tsv]
    :param path_poslist_file: [string]
    :param read_raw_matrix: [bool]
    :return: [annData] annotation data
    """

    # combine paths with filenames
    filtered_matrix_file = path_filtered_files + filenames[0]
    filtered_features_file = path_filtered_files + filenames[1]
    filtered_barcodes_file = path_filtered_files + filenames[2]

    # spatial information file
    spot_locations = pd.read_csv(path_poslist_file, header=None, sep=',')

    # get spatial informations
    spot_locations.rename(columns={0: 'barcode', 1: 'in_tissue', 2: 'array_row', 3: 'array_col',
                                   4: 'pxl_col_in_full_res', 5: 'pxl_row_in_full_res'}, inplace=True)
    spot_locations = spot_locations.sort_values('barcode')

    if read_raw_matrix:
        # combine paths with filenames
        raw_matrix_file = path_raw_files + filenames[0]
        raw_features_file = path_raw_files + filenames[1]
        raw_barcodes_file = path_raw_files + filenames[2]

        # 1. Load RAW data
        raw_adata = sc.read(raw_matrix_file, cache=True)
        raw_adata = raw_adata.transpose()
        # store count matrix in X key of annData
        raw_adata.X = raw_adata.X.toarray()

        raw_file_barcodes = pd.read_csv(raw_barcodes_file, header=None, sep='\t')
        raw_file_features = pd.read_csv(raw_features_file, header=None, sep='\t')

        # # Annotate data
        raw_file_barcodes.rename(columns={0: 'barcode'}, inplace=True)
        raw_file_barcodes.set_index('barcode', inplace=True)
        raw_adata.obs = raw_file_barcodes
        sample_tmp = raw_matrix_file.split(os.sep)[-4]
        raw_adata.obs['sample'] = [sample_tmp] * raw_adata.n_obs
        #   donor = Internal NGI sample indentifier --> have to look up donor (biobank number) and assign it by my own
        raw_adata.obs['project'] = [sample_tmp.split("_")[0]] * raw_adata.n_obs
        #   region = sample number and/or capture area in spatial transcriptomics
        raw_adata.obs['slide'] = [sample_tmp.split("_")[1]] * raw_adata.n_obs

        #   have to additionally include tag_keys for spatial transcriptomics data ..
        raw_file_features.rename(columns={0: 'gene_id', 1: 'gene_name', 2: 'feature_type'}, inplace=True)
        raw_file_features.set_index('gene_name', inplace=True)
        raw_adata.var = raw_file_features
        raw_adata.var_names_make_unique()

        # Binary, indicating if the spot falls inside (1) or outside (0) of tissue.
        raw_adata.obs['in_tissue'] = spot_locations['in_tissue'].values
        # The row coordinate of the spot in the array from 0 to 77. The array has 78 rows.
        raw_adata.obs['array_row'] = spot_locations['array_row'].values
        # The column coordinate of the spot in the array. In order to express the orange crate arrangement of the spots,
        # this column index uses even numbers from 0 to 126 for even rows, and odd numbers from 1 to 127 for odd rows.
        # Notice then that each row (even or odd) has 64 spots.
        raw_adata.obs['array_col'] = spot_locations['array_col'].values
        # The column pixel coordinate of the center of the spot in the full resolution image.
        raw_adata.obs['pxl_col_in_full_res'] = spot_locations['pxl_col_in_full_res'].values
        # The row pixel coordinate of the center of the spot in the full resolution image.
        raw_adata.obs['pxl_row_in_full_res'] = spot_locations['pxl_row_in_full_res'].values
    else:
        # workaround that something can be returend
        raw_adata = []

    # 2. Load FILTERED data
    filtered_adata = sc.read(filtered_matrix_file, cache=True)
    filtered_adata = filtered_adata.transpose()
    filtered_adata.X = filtered_adata.X.toarray()

    filtered_file_barcodes = pd.read_csv(filtered_barcodes_file, header=None, sep='\t')
    filtered_file_features = pd.read_csv(filtered_features_file, header=None, sep='\t')

    # # Annotate data
    filtered_file_barcodes.rename(columns={0: 'barcode'}, inplace=True)
    filtered_file_barcodes.set_index('barcode', inplace=True)
    filtered_adata.obs = filtered_file_barcodes
    sample_tmp = filtered_matrix_file.split(os.sep)[-4]
    filtered_adata.obs['sample'] = [sample_tmp] * filtered_adata.n_obs
    #   donor = Internal NGI sample indentifier --> have to look up donor (biobank number) and assign it by my own
    filtered_adata.obs['project'] = [sample_tmp.split("_")[0]] * filtered_adata.n_obs
    #   region = sample number and/or capture area in spatial transcriptomics todo get slide number from image names
    filtered_adata.obs['slide'] = [sample_tmp.split("_")[1]] * filtered_adata.n_obs
    filtered_adata.obs_names_make_unique()

    # -- get in_tissue index where value == 1
    index_bc_under_tissue = np.where(np.in1d(spot_locations['in_tissue'].values, 1))[0]
    # Binary, indicating if the spot falls inside (1) or outside (0) of tissue.
    filtered_adata.obs['in_tissue'] = spot_locations['in_tissue'].iloc[index_bc_under_tissue].values
    # The row coordinate of the spot in the array from 0 to 77. The array has 78 rows.
    filtered_adata.obs['array_row'] = spot_locations['array_row'].iloc[index_bc_under_tissue].values
    # The column coordinate of the spot in the array. In order to express the orange crate arrangement of the spots,
    # this column index uses even numbers from 0 to 126 for even rows, and odd numbers from 1 to 127 for odd rows.
    # Notice then that each row (even or odd) has 64 spots.
    filtered_adata.obs['array_col'] = spot_locations['array_col'].iloc[index_bc_under_tissue].values
    # The column pixel coordinate of the center of the spot in the full resolution image.
    filtered_adata.obs['pxl_col_in_full_res'] = spot_locations['pxl_col_in_full_res'].iloc[index_bc_under_tissue].values
    # The row pixel coordinate of the center of the spot in the full resolution image.
    filtered_adata.obs['pxl_row_in_full_res'] = spot_locations['pxl_row_in_full_res'].iloc[index_bc_under_tissue].values

    #   have to additionally include tag_keys for spatial transcriptomics data ..
    filtered_file_features.rename(columns={0: 'gene_id', 1: 'gene_name', 2: 'feature_type'}, inplace=True)
    filtered_file_features.set_index('gene_name', inplace=True)
    filtered_adata.var = filtered_file_features
    filtered_adata.var_names_make_unique()

    return raw_adata, filtered_adata


def _scanpy_load_database(path_filtered_files, filenames, path_h5_filtered_matrix, path_h5_raw_matrix, path_spatial,
                          annotation_file_path, sample_id, slide_name, velocyto_path, object_slide, patient,
                          path_raw_files=None, read_raw_matrix=False):
    """
    source: https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/
    Case-study_Mouse-intestinal-epithelium_1906.ipynb

    :param path_raw_files: [string]
        contains path to raw mtx and tsv files created by 10x Genomics Spaceranger
    :param path_filtered_files: [string]
        contains path to filtered mtx and tsv files created by 10x Genomics Spaceranger
    :param filenames: [list]
        contains the names of the mtx and tsv files in the following order:
        [matrix.mtx, features.tsv, barcodes.tsv]
    :param path_h5_filtered_matrix: [string]
    :param path_h5_raw_matrix: [string]
    :param annotation_file_path: [string]
    :param sample_id: [string]
    :param slide_name: [string]
    :param patient: [string]
    :param velocyto_path: [string]
    :param object_slide: [string]
    :param path_spatial: [string]
    :param read_raw_matrix: [bool]
    :return: [annData]
    """

    # combine paths with filenames
    filtered_matrix_file = os.path.join(path_filtered_files, filenames[0])
    filtered_features_file = os.path.join(path_filtered_files, filenames[1])
    filtered_barcodes_file = os.path.join(path_filtered_files, filenames[2])

    if read_raw_matrix:
        # combine paths with filenames
        raw_matrix_file = os.path.join(path_raw_files, filenames[0])
        raw_features_file = os.path.join(path_raw_files, filenames[1])
        raw_barcodes_file = os.path.join(path_raw_files, filenames[2])

        # 1. Load RAW data
        # raw_adata = sc.read_10x_h5(path_h5_matrix, genome=None)
        raw_adata = sc.read(raw_matrix_file, cache=True)
        raw_adata = raw_adata.transpose()
        # store count matrix in X key of annData
        raw_adata.X = raw_adata.X.toarray()

        raw_file_barcodes = pd.read_csv(raw_barcodes_file, header=None, sep='\t')
        raw_file_features = pd.read_csv(raw_features_file, header=None, sep='\t')

        # # Annotate data
        raw_file_barcodes.rename(columns={0: 'barcode'}, inplace=True)
        raw_file_barcodes.set_index('barcode', inplace=True)
        raw_adata.obs = raw_file_barcodes
        sample_tmp = raw_matrix_file.split(os.sep)[-4]
        raw_adata.obs['sample'] = [sample_tmp] * raw_adata.n_obs
        #   donor = Internal NGI sample indentifier --> have to look up donor (biobank number) and assign it by my own
        raw_adata.obs['project'] = [sample_tmp.split("_")[0]] * raw_adata.n_obs
        #   region = sample number and/or capture area in spatial transcriptomics todo get slide number from image names
        # raw_adata.obs['slide'] = [sample_tmp.split("_")[1]] * raw_adata.n_obs
        raw_adata.obs['capture_area'] = [slide_name.split(".")[0]] * raw_adata.n_obs
        # Add object slide serial number
        raw_adata.obs['object_slide'] = [object_slide] * raw_adata.n_obs
        # add patient
        raw_adata.obs['patient'] = [patient] * raw_adata.n_obs
        raw_adata.obs_names_make_unique()

        # include spatial information
        raw_adata = _read_visium(raw_adata, path=[path_h5_raw_matrix, path_spatial])

        # include annotations (no annotations for raw file)
        # raw_adata = _read_excel(raw_adata, path_excel_file=annotation_file_path)

        #   have to additionally include tag_keys for spatial transcriptomics data ..
        raw_file_features.rename(columns={0: 'gene_id', 1: 'gene_name', 2: 'feature_type'}, inplace=True)
        raw_file_features.set_index('gene_name', inplace=True)
        raw_adata.var = raw_file_features
        raw_adata.var_names_make_unique()

        # include spliced, unspliced, Chromosome, End, Start, Strand info (.loom files)
        ldata = scv.read(velocyto_path, cache=True)
        raw_adata = scv.utils.merge(raw_adata, ldata)

        # raw_adata.uns_names_make_unique()
    else:
        # workaround that something can be returned
        raw_adata = []

    # 2. Load FILTERED data
    # filtered_adata = sc.read_10x_h5(path_h5_filtered_matrix, genome=None)

    # if 'P21093_21L008958' in filtered_matrix_file:
    #     print('stop')
    filtered_adata = sc.read(filtered_matrix_file, cache=True)
    filtered_adata = filtered_adata.transpose()
    filtered_adata.X = filtered_adata.X.toarray()

    filtered_file_barcodes = pd.read_csv(filtered_barcodes_file, header=None, sep='\t')
    filtered_file_features = pd.read_csv(filtered_features_file, header=None, sep='\t')

    # # Annotate data
    filtered_file_barcodes.rename(columns={0: 'barcode'}, inplace=True)
    filtered_file_barcodes.set_index('barcode', inplace=True)
    filtered_adata.obs = filtered_file_barcodes
    sample_tmp = filtered_matrix_file.split(os.sep)[-4]
    filtered_adata.obs['sample'] = [sample_tmp] * filtered_adata.n_obs
    #   donor = Internal NGI sample indentifier --> have to look up donor (biobank number) and assign it by my own
    filtered_adata.obs['project'] = [sample_tmp.split("_")[0]] * filtered_adata.n_obs
    #   region = sample number and/or capture area in spatial transcriptomics
    # filtered_adata.obs['slide'] = [sample_tmp.split("_")[1]] * filtered_adata.n_obs
    filtered_adata.obs['capture_area'] = [slide_name.split(".")[0]] * filtered_adata.n_obs
    # Add object slide serial number
    filtered_adata.obs['object_slide'] = [object_slide] * filtered_adata.n_obs
    # add patient
    filtered_adata.obs['patient'] = [patient] * filtered_adata.n_obs
    filtered_adata.obs_names_make_unique()

    # include spatial information
    filtered_adata = _read_visium(adata=filtered_adata, path=[path_h5_filtered_matrix, path_spatial],
                                  library_id=sample_id)

    # include annotations (no info for project P16357)
    filtered_adata = _read_excel(filtered_adata, path_excel_file=annotation_file_path)

    #   have to additionally include tag_keys for spatial transcriptomics data ..
    filtered_file_features.rename(columns={0: 'gene_id', 1: 'gene_name', 2: 'feature_type'}, inplace=True)
    filtered_file_features.set_index('gene_name', inplace=True)
    filtered_adata.var = filtered_file_features
    filtered_adata.var_names_make_unique()

    # include spliced, unspliced, Chromosome, End, Start, Strand info (.loom files)
    ldata = scv.read(velocyto_path, cache=True)
    filtered_adata = scv.utils.merge(filtered_adata, ldata)

    return raw_adata, filtered_adata


def main(filename, read_raw_matrix=False, spatial_concat=True):
    """
    Load samples with count matrix, features, barcodes, spatial information,
    Additionally add the information of manual annotations and velocyto

    :param filename: [string]
        path to be configuration file to read out data sets (type string)
    :param read_raw_matrix: [bool]
        only true if you really want to load the giant unfiltered count matrix
        if raw matrix shall be loaded for further analysis
        -> can be used to study diffusion outside of tissue or later diffusion inside the tissue
    :param spatial_concat: [bool] if more than 1 sample and samples contain spatial information --> True
        if samples shall be saved together in a annData object or in a list
    :return: [annData]
        raw and filtered read out matrices and spatial information and annotation dataset
    """
    # check whether to do a single sample load or if you have more than one sample to load
    # input_path = os.path.join(os.environ['PYTHONPATH'].split(os.pathsep)[0], 'Input', 'config_files', filename)
    config_paths = ht.load_sample_config_file(filename=filename, file_type="csv")

    absolute_path = os.path.join(os.sep, config_paths[0][0], config_paths[2][0])
    # matrix_file_end, features_file_end, barcode_file_end
    feature_bc_matrix_string = np.array([config_paths[9][0], config_paths[8][0], config_paths[7][0]])

    # # Parse Filenames
    project = config_paths[2][0]  # if more examples then use sample_strings.pop(0)
    sample_id = config_paths[3][0]
    excel_sheet_id = config_paths[20][0]  # if more examples then use config_paths[20].pop(0)
    loom_id = config_paths[22][0]  # contains the information about spliced and unspliced genes

    # save sample ids in list
    list_sample_ids = [sample_id]

    if read_raw_matrix:
        print("\nRaw/Unfiltered feature-barcode matrix contains every barcode from fixed list of known-good barcode "
              "sequences. This includes background and cell associated barcodes")

        # path to raw files ending with .mtx and .tsv (type string)
        raw_feature_bc_matrix_path = os.path.join(absolute_path, config_paths[1][0], project + "_" + sample_id,
                                                  config_paths[4][0], config_paths[5][0])
        print("Filtered feature-barcode matrix contains only cells associated barcodes")
        # path h5
        filtered_bc_matrix_h5_path = os.path.join(absolute_path, config_paths[1][0], project + "_" + sample_id,
                                                  config_paths[4][0], config_paths[11][0])
        raw_bc_matrix_h5_path = os.path.join(absolute_path, config_paths[1][0], project + "_" + sample_id,
                                             config_paths[4][0], config_paths[10][0])

        # path to filtered files ending with .mtx and .tsv (type string)
        filtered_feature_bc_matrix_path = os.path.join(absolute_path, config_paths[1][0], project + "_" + sample_id,
                                                       config_paths[4][0], config_paths[6][0])

        # path to spatial information
        spatial_path = os.path.join(absolute_path, config_paths[1][0], project + "_" + sample_id, config_paths[12][0])

        # path to annotation files
        annotation_file_path = os.path.join(absolute_path, config_paths[19][0], excel_sheet_id)

        # path to velocyto information
        velo_path = os.path.join(config_paths[0][0], config_paths[2][0], config_paths[21][0],
                                 project + "_" + sample_id, loom_id)

        # # Annotate data
        print("\n-------- Start: Read out values --------")
        # Two options to read in feature_ids, gene_names, feature_types, barcodes, count_matrix_data
        # 1. Malte Luecken using Scanpy from TheisLab; read out mtx and tsv files
        raw_annot_data, filtered_annot_data = _scanpy_load_database(
            path_filtered_files=filtered_feature_bc_matrix_path,
            path_raw_files=raw_feature_bc_matrix_path,
            path_h5_filtered_matrix=filtered_bc_matrix_h5_path,
            path_h5_raw_matrix=raw_bc_matrix_h5_path,
            filenames=feature_bc_matrix_string,
            path_spatial=spatial_path,
            annotation_file_path=annotation_file_path,
            read_raw_matrix=read_raw_matrix,
            sample_id=project + "_" + sample_id,
            slide_name=config_paths[16][0],
            velocyto_path=velo_path,
            object_slide=config_paths[23][0],
            patient=config_paths[24][0])

        # initialize list of adatas
        list_filtered_annot_data = [filtered_annot_data]
        list_raw_annot_adata = [raw_annot_data]

        for c_sample in tqdm(range(len(config_paths[0][1:])), desc='Loading samples'):
            c_sample += 1
            # path to h5 and matrices
            path_matrix = os.path.join(os.sep, config_paths[0][c_sample], config_paths[2][c_sample],
                                       config_paths[1][c_sample], config_paths[3][c_sample], config_paths[4][c_sample])

            # path h5
            filtered_bc_matrix_h5_path = os.path.join(path_matrix, config_paths[11][c_sample])
            raw_bc_matrix_h5_path = os.path.join(path_matrix, config_paths[10][c_sample])

            # matrix
            raw_feature_bc_matrix_path = os.path.join(path_matrix, config_paths[5][c_sample])
            filtered_feature_bc_matrix_path = os.path.join(path_matrix, config_paths[6][c_sample])

            #  path to spatial information
            spatial_path = os.path.join(os.sep, config_paths[0][c_sample], config_paths[2][c_sample],
                                        config_paths[12][c_sample])

            # path to annotation files
            annotation_file_path = os.path.join(os.sep, config_paths[0][c_sample], config_paths[2][c_sample],
                                                config_paths[19][c_sample], config_paths[20][c_sample])

            # path to velocyto information
            velo_path = os.path.join(os.sep, config_paths[0][c_sample], config_paths[2][c_sample],
                                     config_paths[21][c_sample],
                                     config_paths[2][c_sample] + "_" + config_paths[3][c_sample],
                                     config_paths[23][c_sample])

            # # Load spatial information and count matrix
            raw_adata_tmp, filtered_adata_tmp = _scanpy_load_database(
                path_filtered_files=filtered_feature_bc_matrix_path,
                path_raw_files=raw_feature_bc_matrix_path,
                path_h5_filtered_matrix=filtered_bc_matrix_h5_path,
                path_h5_raw_matrix=raw_bc_matrix_h5_path,
                path_spatial=spatial_path,
                annotation_file_path=annotation_file_path,
                read_raw_matrix=read_raw_matrix,
                sample_id=config_paths[2][c_sample] + "_" + config_paths[3][c_sample],
                slide_name=config_paths[16][c_sample], filenames=feature_bc_matrix_string,
                velocyto_path=velo_path, object_slide=config_paths[23][c_sample],
                patient=config_paths[24][c_sample])

            # # Concatenate data sets (also do this if you have more than one donor!)
            if spatial_concat:
                #  FILTERED
                list_filtered_annot_data.append(filtered_adata_tmp)

                #   RAW
                list_raw_annot_adata.append(raw_adata_tmp)
            else:
                #   RAW
                raw_annot_data = raw_annot_data.concatenate(raw_adata_tmp, batch_key='sample_id')
                raw_annot_data.obs.drop(columns=['sample_id'], inplace=True)
                raw_annot_data.obs_names = [c.split("-")[0] for c in raw_annot_data.obs_names]
                raw_annot_data.obs_names_make_unique()

                #  FILTERED
                filtered_annot_data = filtered_annot_data.concatenate(filtered_adata_tmp, batch_key='sample_id')
                filtered_annot_data.obs.drop(columns=['sample_id'], inplace=True)
                filtered_annot_data.obs_names = [c.split("-")[0] for c in filtered_annot_data.obs_names]
                filtered_annot_data.obs_names_make_unique()

            # save sample ids in list
            list_sample_ids.append(config_paths[3][c_sample])

        # workaround because concatenate cannot handle concat and vanilla adata
        if spatial_concat:
            filtered_annot_data = ht.concatenate_adatas(list_filtered_annot_data)
            # make barcode names unique
            filtered_annot_data.obs_names = [c.split("-")[0] for c in filtered_annot_data.obs_names]

            # RAW
            raw_annot_data = ht.concatenate_adatas(list_raw_annot_adata)
            # make barcode names unique
            raw_annot_data.obs_names = [c.split("-")[0] for c in raw_annot_data.obs_names]
    else:
        print("Filtered feature-barcode matrix contains only spots associated barcodes under tissue")
        # path to h5 and matrices
        path_matrix = os.path.join(os.sep, config_paths[0][0], config_paths[2][0], config_paths[1][0],
                                   config_paths[2][0] + "_" + config_paths[3][0], config_paths[4][0])

        # path h5
        filtered_bc_matrix_h5_path = os.path.join(path_matrix, config_paths[11][0])
        raw_bc_matrix_h5_path = os.path.join(path_matrix, config_paths[10][0])

        # path to filtered files ending with .mtx and .tsv (type string)
        filtered_feature_bc_matrix_path = os.path.join(path_matrix, config_paths[6][0])

        # path to spatial information
        spatial_path = os.path.join(os.sep, config_paths[0][0], config_paths[2][0], config_paths[1][0],
                                    config_paths[2][0] + "_" + config_paths[3][0], config_paths[12][0])

        # path to annotation files
        annotation_file_path = os.path.join(os.sep, config_paths[0][0], config_paths[2][0], config_paths[19][0],
                                            excel_sheet_id)

        # path to velocyto information
        velo_path = os.path.join(os.sep, config_paths[0][0], config_paths[2][0], config_paths[21][0],
                                 project + "_" + sample_id, loom_id)

        # # Annotate data
        print("\n-------- Start: Read out values --------")
        # Two options to read in feature_ids, gene_names, feature_types, barcodes, count_matrix_data
        # 1. Malte Luecken using Scanpy from TheisLab; read out mtx and tsv files
        raw_annot_data, filtered_annot_data = _scanpy_load_database(
            path_filtered_files=filtered_feature_bc_matrix_path,
            filenames=feature_bc_matrix_string,
            path_spatial=spatial_path,
            path_h5_filtered_matrix=filtered_bc_matrix_h5_path,
            path_h5_raw_matrix=raw_bc_matrix_h5_path,
            annotation_file_path=annotation_file_path,
            sample_id=config_paths[2][0] + "_" + sample_id,
            slide_name=config_paths[16][0],
            velocyto_path=velo_path, object_slide=config_paths[23][0],
            patient=config_paths[24][0])

        # initialize list of adatas
        list_filtered_annot_data = [filtered_annot_data]
        list_raw_annot_adata = [raw_annot_data]

        # # Loop to load all data sets
        for c_sample in tqdm(range(len(config_paths[0][1:])), desc='Loading samples'):
            c_sample += 1
            print("_".join([config_paths[2][c_sample], config_paths[3][c_sample]]))
            # path to h5 and matrices
            path_matrix = os.path.join(os.sep, config_paths[0][c_sample], config_paths[2][c_sample],
                                       config_paths[1][c_sample],
                                       config_paths[2][c_sample] + "_" + config_paths[3][c_sample],
                                       config_paths[4][c_sample])

            # path h5
            filtered_bc_matrix_h5_path = os.path.join(path_matrix, config_paths[11][c_sample])
            raw_bc_matrix_h5_path = os.path.join(path_matrix, config_paths[10][c_sample])

            filtered_feature_bc_matrix_path = os.path.join(path_matrix, config_paths[6][c_sample])

            #  path to spatial information
            spatial_path = os.path.join(os.sep, config_paths[0][c_sample], config_paths[2][c_sample],
                                        config_paths[1][c_sample],
                                        config_paths[2][c_sample] + "_" + config_paths[3][c_sample],
                                        config_paths[12][c_sample])

            # path to annotation files
            annotation_file_path = os.path.join(os.sep, config_paths[0][c_sample], config_paths[2][c_sample],
                                                config_paths[19][c_sample], config_paths[20][c_sample])

            # path to velocyto information
            velo_path = os.path.join(os.sep, config_paths[0][c_sample], config_paths[2][c_sample],
                                     config_paths[21][c_sample],
                                     config_paths[2][c_sample] + "_" + config_paths[3][c_sample],
                                     config_paths[22][c_sample])

            # # Load spatial information and count matrix
            raw_adata_tmp, filtered_adata_tmp = _scanpy_load_database(
                path_filtered_files=filtered_feature_bc_matrix_path,
                filenames=feature_bc_matrix_string,
                path_spatial=spatial_path,
                path_h5_filtered_matrix=filtered_bc_matrix_h5_path,
                path_h5_raw_matrix=raw_bc_matrix_h5_path,
                annotation_file_path=annotation_file_path,
                sample_id=config_paths[2][c_sample] + "_" + config_paths[3][c_sample],
                slide_name=config_paths[16][c_sample],
                velocyto_path=velo_path, object_slide=config_paths[23][c_sample],
                patient=config_paths[24][c_sample])

            # # Concatenate data sets (also do this if you have more than one donor!)
            if spatial_concat:
                list_filtered_annot_data.append(filtered_adata_tmp)
                list_raw_annot_adata.append(raw_adata_tmp)
            else:
                filtered_annot_data = filtered_annot_data.concatenate(filtered_adata_tmp, batch_key='sample_id')
                # filtered_annot_data.var['gene_id'] = filtered_annot_data.var['gene_id-1']
                # filtered_annot_data.var.drop(columns=['gene_id-1', 'gene_id-0'], inplace=True)
                filtered_annot_data.obs.drop(columns=['sample_id'], inplace=True)
                filtered_annot_data.obs_names = [c.split("-")[0] for c in filtered_annot_data.obs_names]
                # make filtered_annot_data.obs_names only unique if you have scRNA-seq data set

            # save sample ids in list
            list_sample_ids.append(config_paths[3][c_sample])

        # workaround because concatenate cannot handle concat and vanilla adata
        if spatial_concat:
            filtered_annot_data = ht.concatenate_adatas(list_filtered_annot_data)
            # make barcode names unique
            filtered_annot_data.obs_names = [c.split("-")[0] for c in filtered_annot_data.obs_names]

    return raw_annot_data, filtered_annot_data, config_paths, list_sample_ids
