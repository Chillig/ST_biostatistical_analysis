import json
import numpy as np
import pandas as pd
import scanpy as sc


def _scanpy_load_annotate_tsv_mtx_files(path_filtered_files, filenames, path_raw_files=None, read_raw_matrix=False):
    """
    source: https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/
    Case-study_Mouse-intestinal-epithelium_1906.ipynb

    :param path_raw_files: contains path to raw mtx and tsv files created by 10x Genomics Spaceranger
    :param path_filtered_files: contains path to filtered mtx and tsv files created by 10x Genomics Spaceranger
    :param filenames: contains the names of the mtx and tsv files in the following order:
    [matrix.mtx, features.tsv, barcodes.tsv]
    :param read_raw_matrix:
    :return: annotation data
    """
    filtered_matrix_file = path_filtered_files + filenames[0]
    filtered_features_file = path_filtered_files + filenames[1]
    filtered_barcodes_file = path_filtered_files + filenames[2]

    if read_raw_matrix:
        # combine paths with filenames
        raw_matrix_file = path_raw_files + filenames[0]
        raw_features_file = path_raw_files + filenames[1]
        raw_barcodes_file = path_raw_files + filenames[2]

        # 1. Load RAW data
        raw_adata = sc.read(raw_matrix_file)
        raw_adata = raw_adata.transpose()
        # store count matrix in X key of annData
        raw_adata.X = raw_adata.X.toarray()

        raw_file_barcodes = pd.read_csv(raw_barcodes_file, header=None, sep='\t')
        raw_file_features = pd.read_csv(raw_features_file, header=None, sep='\t')

        # # Annotate data
        raw_file_barcodes.rename(columns={0: 'barcode'}, inplace=True)
        raw_file_barcodes.set_index('barcode', inplace=True)
        raw_adata.obs = raw_file_barcodes
        sample_tmp = path_raw_files.split("/")[-4]
        raw_adata.obs['sample'] = [sample_tmp] * raw_adata.n_obs
        #   donor = Internal NGI sample indentifier --> have to look up donor (biobank number) and assign it by my own
        raw_adata.obs['donor'] = [sample_tmp[:-1]] * raw_adata.n_obs
        #   region = sample number and/or capture area in spatial transcriptomics
        raw_adata.obs['region'] = [sample_tmp[-1]] * raw_adata.n_obs

        #   have to additionally include tag_keys for spatial transcriptomics data ..
        raw_file_features.rename(columns={0: 'gene_id', 1: 'gene_name', 2: 'feature_type'}, inplace=True)
        raw_file_features.set_index('gene_name', inplace=True)
        raw_adata.var = raw_file_features
        raw_adata.var_names_make_unique()
    else:
        # workaround that something can be returend
        raw_adata = []

    # 2. Load FILTERED data
    filtered_adata = sc.read(filtered_matrix_file)
    filtered_adata = filtered_adata.transpose()
    filtered_adata.X = filtered_adata.X.toarray()

    filtered_file_barcodes = pd.read_csv(filtered_barcodes_file, header=None, sep='\t')
    filtered_file_features = pd.read_csv(filtered_features_file, header=None, sep='\t')

    # # Annotate data
    filtered_file_barcodes.rename(columns={0: 'barcode'}, inplace=True)
    filtered_file_barcodes.set_index('barcode', inplace=True)
    filtered_adata.obs = filtered_file_barcodes
    sample_tmp = path_filtered_files.split("/")[-4]
    filtered_adata.obs['sample'] = [sample_tmp] * filtered_adata.n_obs
    #   donor = Internal NGI sample indentifier --> have to look up donor (biobank number) and assign it by my own
    filtered_adata.obs['donor'] = [sample_tmp[:-1]] * filtered_adata.n_obs
    #   region = sample number and/or capture area in spatial transcriptomics
    filtered_adata.obs['region'] = [sample_tmp[-1]] * filtered_adata.n_obs
    filtered_adata.obs_names_make_unique()

    #   have to additionally include tag_keys for spatial transcriptomics data ..
    filtered_file_features.rename(columns={0: 'gene_id', 1: 'gene_name', 2: 'feature_type'}, inplace=True)
    filtered_file_features.set_index('gene_name', inplace=True)
    filtered_adata.var = filtered_file_features
    filtered_adata.var_names_make_unique()

    return raw_adata, filtered_adata


def load_config_file(json_filename):
    """
    Load config file with paths to counts matrix and additonal information

    :param json_filename:
    :return:
    """

    with open(json_filename) as json_file:
        configs = json.load(json_file)

        # No need to save configs in different variables; it is just for visualization and easier to read out later
        file_path = configs['file_path']
        sample_strings = configs['sample_strings']
        sample_id_strings = configs['sample_id_strings']

        # output matrix
        output_type_matrix = configs['output_type_matrix']
        raw_feature_bc_matrix_folder = configs['raw_feature_bc_matrix_folder']
        filtered_feature_bc_matrix_folder = configs['filtered_feature_bc_matrix_folder']

        barcode_file_end = configs['barcode_file_end']
        features_file_end = configs['features_file_end']
        matrix_file_end = configs['matrix_file_end']
        raw_hdf5_file_end = configs['raw_hdf5_file_end']
        filtered_hdf5_file_end = configs['filtered_hdf5_file_end']

    return \
        file_path, sample_strings, sample_id_strings, output_type_matrix, raw_feature_bc_matrix_folder, \
        filtered_feature_bc_matrix_folder, barcode_file_end, features_file_end, matrix_file_end, raw_hdf5_file_end,\
        filtered_hdf5_file_end


def main(json_filename, read_raw_matrix=False):
    """

    :param json_filename: path to be configuration file to read out data sets (type string)
    :param read_raw_matrix: only true if you really want to load the giant unfiltered count matrix
    :return: raw and filtered read out matrices and annotation dataset
    """

    config_paths = load_config_file(json_filename=json_filename)

    # matrix_file_end, features_file_end, barcode_file_end
    feature_bc_matrix_string = np.array([config_paths[8], config_paths[7], config_paths[6]])

    # # Parse Filenames
    sample = config_paths[1]  # if more examples then use sample_strings.pop(0)
    sample_id = config_paths[2].pop(0)

    if read_raw_matrix:
        print("\nRaw/Unfiltered feature-barcode matrix contains every barcode from fixed list of known-good barcode "
              "sequences. This includes background and cell associated barcodes")

        # path to raw files ending with .mtx and .tsv (type string)
        raw_feature_bc_matrix_path = config_paths[0] + sample + sample_id + config_paths[3] + config_paths[4]
        print("Filtered feature-barcode matrix contains only cells associated barcodes")
        # path to filtered files ending with .mtx and .tsv (type string)
        filtered_feature_bc_matrix_path = config_paths[0] + sample + sample_id + config_paths[3] + config_paths[5]

        # # Annotate data
        print("\n-------- Start: Read out values --------")
        # Two options to read in feature_ids, gene_names, feature_types, barcodes, count_matrix_data
        # 1. Malte Luecken using Scanpy from TheisLab; read out mtx and tsv files
        raw_annot_data, filtered_annot_data = _scanpy_load_annotate_tsv_mtx_files(raw_feature_bc_matrix_path,
                                                                                  filtered_feature_bc_matrix_path,
                                                                                  feature_bc_matrix_string,
                                                                                  read_raw_matrix=read_raw_matrix)

        # # Loop to load all data sets
        for i_counter in range(len(config_paths[2])):
            print("         ... reading out ...")
            raw_feature_bc_matrix_path = \
                config_paths[0] + sample + config_paths[2][i_counter] + config_paths[3] + config_paths[4]
            filtered_feature_bc_matrix_path = \
                config_paths[0] + sample + config_paths[2][i_counter] + config_paths[3] + config_paths[5]

            # # Load count matrix, features and observables
            raw_adata_tmp, filtered_adata_tmp = _scanpy_load_annotate_tsv_mtx_files(raw_feature_bc_matrix_path,
                                                                                    filtered_feature_bc_matrix_path,
                                                                                    feature_bc_matrix_string,
                                                                                    read_raw_matrix=read_raw_matrix)

            # # Concatenate data sets (also do this if you have more than one donor!)
            #   RAW
            raw_annot_data = raw_annot_data.concatenate(raw_adata_tmp, batch_key='sample_id')
            # raw_annot_data.var['gene_id'] = raw_annot_data.var['gene_id-1']
            # raw_annot_data.var.drop(columns=['gene_id-1', 'gene_id-0'], inplace=True)
            raw_annot_data.obs.drop(columns=['sample_id'], inplace=True)
            raw_annot_data.obs_names = [c.split("-")[0] for c in raw_annot_data.obs_names]
            raw_annot_data.obs_names_make_unique(join='_')
            # join '_' does not work here: barcodes are the same in each capture area!!
            # raw_annot_data.obs_names_make_unique(join='_')

            #  FILTERED
            filtered_annot_data = filtered_annot_data.concatenate(filtered_adata_tmp, batch_key='sample_id')
            # filtered_annot_data.var['gene_id'] = filtered_annot_data.var['gene_id-1']
            # filtered_annot_data.var.drop(columns=['gene_id-1', 'gene_id-0'], inplace=True)
            filtered_annot_data.obs.drop(columns=['sample_id'], inplace=True)
            filtered_annot_data.obs_names = [c.split("-")[0] for c in filtered_annot_data.obs_names]
            filtered_annot_data.obs_names_make_unique(join='_')
            # join '_' does not work here: barcodes are the same in each capture area!!
            # filtered_annot_data.obs_names_make_unique(join='_')

    else:
        print("Filtered feature-barcode matrix contains only cells associated barcodes")
        # path to filtered files ending with .mtx and .tsv (type string)
        filtered_feature_bc_matrix_path = config_paths[0] + sample + sample_id + config_paths[3] + config_paths[5]

        # # Annotate data
        print("\n-------- Start: Read out values --------")
        # Two options to read in feature_ids, gene_names, feature_types, barcodes, count_matrix_data
        # 1. Malte Luecken using Scanpy from TheisLab; read out mtx and tsv files
        _, filtered_annot_data = _scanpy_load_annotate_tsv_mtx_files(filtered_feature_bc_matrix_path,
                                                                     feature_bc_matrix_string)

        # # Loop to load all data sets
        for i_counter in range(len(config_paths[2])):
            print("         ... reading out ...")
            filtered_feature_bc_matrix_path = \
                config_paths[0] + sample + config_paths[2][i_counter] + config_paths[3] + config_paths[5]

            # # Load count matrix, features and observables
            _, filtered_adata_tmp = _scanpy_load_annotate_tsv_mtx_files(filtered_feature_bc_matrix_path,
                                                                        feature_bc_matrix_string)

            #  FILTERED
            filtered_annot_data = filtered_annot_data.concatenate(filtered_adata_tmp, batch_key='sample_id')
            filtered_annot_data.obs.drop(columns=['sample_id'], inplace=True)
            filtered_annot_data.obs_names = [c.split("-")[0] for c in filtered_annot_data.obs_names]
            filtered_annot_data.obs_names_make_unique()

        # Workaround to ensure that something can be returned
        raw_annot_data = []

    # ---- Side notes to know but can also be looked up in the summary created by 10x Genomics Spaceranger --- #
    unique_sample = np.unique(filtered_annot_data.obs["sample"])
    num_cells_previous_sample = 0
    for sample_name in unique_sample:
        print("\nSample {}: ".format(sample_name))
        print("\nSide notes of {} ".format(filtered_annot_data.obs['sample'].values[1]))
        number_cells = len(np.where(filtered_annot_data.obs['sample'].values == sample_name)[0])

        print("No. cells: ", number_cells)
        # Count number of expressed genes (count one gene over all spots)
        counts_gene = filtered_annot_data[num_cells_previous_sample:number_cells].X.sum(0)
        counts_gene_sorted = np.sort(counts_gene)
        print("Total No. genes detected: ", np.count_nonzero(counts_gene_sorted))

        # Calculate median genes per spot
        copy_s1 = filtered_annot_data[:number_cells].X.copy()
        mask = copy_s1 > 0
        zero_array = np.zeros_like(copy_s1)
        # count numbers of True == numbers of gene overall spots
        zero_array[mask] = 1
        median_genes_per_spot = np.median(zero_array.sum(1))
        median_umi_counts_per_spot = np.median(copy_s1.sum(1))
        print("Median genes: ", median_genes_per_spot)
        print("Total No. of UMI Counts: ", sum(copy_s1.sum(1)))
        print("Median UMI Counts: ", median_umi_counts_per_spot)
        num_cells_previous_sample = number_cells

    # All samples
    # Get barcodes of each sample and therefore the number of cells for each sample
    model = filtered_annot_data.obs[['sample'] + []]
    batch_info = model.groupby('sample').groups.values()
    n_batches = np.array([len(v) for v in batch_info])
    print("\n Sorted No. genes for each sample: ", n_batches)

    print("\n")
    # ---                                            End side notes                                         --- #

    # Second option: load from hdf5 files

    return raw_annot_data, filtered_annot_data, config_paths
