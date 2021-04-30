# Python version '3.7.7 (default, Mar 23 2020, 17:31:31) \n[Clang 4.0.1 (tags/RELEASE_401/final)]'
# Platform: Darwin-19.4.0-x86_64-i386-64bit
# Running under: MacOS 10.15.4 Catalina

from datetime import date
import os
import glob
import scanpy as sc
import numpy as np
import pandas as pd

import scripts.pre_processing.main_preprocessing as prep
import Source.Approach_III.main_tasks as dge_approach
import Source.Downstream_analysis.main_downstream_analysis as dsa
import Source.Spatial_velocytos.spatial_velocyto as spatial_velo
from scripts.utils import loading_matrices

# initialization variables
# get todays date
adata_save_date = date.today()  # "2020-10-06"
# todo save all in one folder together with adata object ..
# Path to Database
database_path = \
    'Volumes' + os.sep + 'Samsung_T5' + os.sep + 'IGSSE-TUM__datasets' + os.sep + \
    '2020__Kilian_Eyerich__Spatial_transcriptomics_skin_single-cell_RNA-seq' + os.sep + 'Database'
# decisions analysis type
analysis_sinlge_sample = False
# decision processing steps
pre_processing = True
dge_analysis = False
downstream_analyis = False
spatial_velocyto_analysis = False
# decisions pre-processing:
read_raw_matrix = False
# combined analysis on all samples if from different data sets and samples splitted into batches
spatial_concat = True
apply_mt_threshold = True
norm_type = 'scanpy'
exclude_highly_expressed = False
apply_remove_cc_effect = False
bc_algorithm = 'scanorama_bc'
# ==>  apply first data integration then batch correction!
apply_data_integration_bc = False
cc_genes_file = 'Input' + os.sep + 'cell_cycle_files' + os.sep + 'Macosko_cell_cycle_genes_2015.csv'

# # Data files
config_filename = 'P15509_P16357'
project_name = ['P15509', 'P16357']
# ids for project number P16357
ids = np.arange(1001, 1061, 1)
# all ids from project P15509 and P16357
sample_id = ['1001', '1002', '1003', '1004'] + list(ids.astype(str))

# # DGE Analysis decisions
experiment = "normed_by_spots_pseudo-bulk"
save_supadata = False


def concatenate_list_data(string_list):
    """

    :param string_list: [list]
    :return:
    """
    result = ''
    for element in string_list:
        result += "_" + str(element)
    return result


def save_gene_symbol_ensembl_map(adata, save_folder):
    """
    Create .csv file to map gene symbols to ENSEMBL gene name because no unqiue mapping is possible afterwards
    - Gene definition in ENSEMBL is locus-based because it is associated with the reference genome: here GRCH38
    - Gene definition of HGNC is defined by its contribution to a phenotype or function
    => no unique mapping of gene symbol to ENSEMBL gene

    :param adata: [annData]
    :param save_folder: [string]
    :return:
    """

    df_gene = pd.DataFrame({"symbol": adata.var.index, "ensembl": adata.var["gene_id"].values})
    df_gene.to_csv(os.path.join(save_folder, 'map_gene_symbol_ensembl.csv'))


def save_count_matrix_to_csv(adata, save_folder):
    """

    :param adata: [annData]
    :param save_folder: [string]
    :return:
    """
    observables = list(adata.obs_keys())
    # 1. Get all manual annotations by extracting all keys with upper case characters
    annotations = [char for char in observables if any(c.isupper() for c in char)]

    if "ANNOTATOR" in annotations:
        annotations.remove("ANNOTATOR")

    # 3. Remove cell cycle annotations from tissue_cell_labels list
    for cc in ["G1_score", "G2M_score", "S_score", "M_score"]:
        try:
            annotations.remove(cc)
        except ValueError:
            continue

    samples = np.unique(adata.obs['sample'])

    # save count matrix for each sample: gene x sample to load it in R
    for obs in annotations:
        df = pd.DataFrame()
        for s_id in samples:
            if 'counts' in adata.layers.keys():
                count_matrix = adata.layers['counts'][np.where((adata.obs['sample'] == s_id) &
                                                               (adata.obs[obs] == 1))[0]]
            else:
                count_matrix = adata.X[np.where((adata.obs['sample'] == s_id) & (adata.obs[obs] == 1))[0]]

            gene_expr = count_matrix.sum(axis=0)
            df[s_id] = gene_expr

            # save count matrix barcodes x genes for each sample
            barcodes = adata.obs.index[np.where((adata.obs['sample'] == s_id) & (adata.obs[obs] == 1))[0]]
            df_cm = pd.DataFrame(count_matrix, index=barcodes, columns=adata.var.index)
            df_cm.to_csv(os.path.join(save_folder, "_".join(["_".join([s_id]),  "_".join([obs]),
                                                             'splitted_adata_R.csv'])), sep='\t')
        df = df.set_index(adata.var.index)
        df.to_csv(os.path.join(save_folder, "_".join(['joined', "_".join([obs]), 'adata_R.csv'])), sep='\t')


def dataset_processing(adata, filelist, results_save_folder, adata_filename, filename, single_sample,
                       preprosessing=True, dge=True, downstream=True, spatial_velocyto=True):
    """
    Decide if pre-processing needs to be applied and which kind of analysis steps should be down
    Possible computing steps on annData onject are:
    1. Pre-processing
    2. DGE, Go-term and Pathway Enrichment Analysis
    3. Downstream analysis

    :param adata: [annData]
    :param filelist: [string]
    :param results_save_folder: [string]
    :param adata_filename: [string]
    :param filename: [string]
    :param single_sample: [bool]
    :param preprosessing: [bool]
    :param dge: [bool]
    :param downstream: [bool]
    :param spatial_velocyto: [bool]
    :return:
    """

    if preprosessing:
        print("#   --  >Pre-processing<  --   #")
        adata, adata_pathfile = prep.main(adata=adata, save_folder=results_save_folder,
                                          adata_path_filenames=adata_filename, cc_genes_file=cc_genes_file,
                                          bc_algorithm=bc_algorithm, norm_type=norm_type,
                                          apply_mt_threshold=apply_mt_threshold,
                                          apply_remove_cc_effect=apply_remove_cc_effect,
                                          read_raw_matrix=read_raw_matrix, spatial_concat=spatial_concat,
                                          apply_data_integration_bc=apply_data_integration_bc,
                                          exclude_highly_expressed=exclude_highly_expressed)

    print("Initialise Analysis Steps")
    if dge:
        adata_file = [f for f in filelist if '{}_QC.h5'.format(adata_filename) in f]
        adata_dge = sc.read(adata_file[0])
        print("\n#   --  >DGE, Go-term and Pathway Enrichment Analysis<  --   #")
        dge_approach.main(adata=adata_dge, experiment=experiment, save_supadata=save_supadata,
                          save_folder=results_save_folder, project_id=filename, raw=read_raw_matrix)

    if downstream:
        print("\n#   --  >Downstream Analysis<  --   #")
        dsa.main(adata=adata, save_folder=results_save_folder, file_name=filename)

    if spatial_velocyto:
        print("\n#   --  >Spatial Velocyto Analysis<  --   #")
        spatial_velo.main(adata=adata, save_folder=results_save_folder, single_adata=single_sample)


def main():
    """
    Load annData object(s) with variables, observables, unstructured data ...
    Pre-processing phase: QC + Filtering out biological and technical variations
    Downstream analysis: Clustering, Cell type annotation, DGE, ...


    Set seed for reproducible research if not using scanpy (scanpy does this automatically)
    (stochastic processes change their outcome if no seed is set at beginning random.seed(123))

    """
    if adata_save_date == "today":
        cw_date = date.today()
    else:
        cw_date = adata_save_date

    res_date = date.today()

    current_working_directory = os.environ['PYTHONPATH'].split(os.pathsep)[0]

    # create adata saving place
    saving_adata_dir = os.path.join(current_working_directory, os.path.join('adata_storage', str(cw_date)))
    os.makedirs(saving_adata_dir, exist_ok=True)

    # create plot and files saving place for adata object
    if read_raw_matrix:
        results_save_folder = os.path.join(current_working_directory, os.path.join('Output', 'results', str(res_date),
                                                                                   'raw'))
    else:
        results_save_folder = os.path.join(current_working_directory, os.path.join('Output', 'results', str(res_date),
                                                                                   'filtered'))
    os.makedirs(results_save_folder, exist_ok=True)

    if analysis_sinlge_sample:
        # Do the analysis for each sample separately from one project
        for ind, id_number in enumerate(sample_id):
            # set save_folder path for results
            save_folder_ss_analysis = os.path.join(results_save_folder, config_filename)
            os.makedirs(save_folder_ss_analysis, exist_ok=True)

            # get files in adata storage path if there are non than filelist will be empty []
            filelist = glob.glob(os.path.join(saving_adata_dir,
                                              '*_adata_{}_{}*.h5'.format(project_name, sample_id[ind])))

            # check if files exist
            files_exist = [f for f in filelist if '_pp' in f]
            # files_exist = [f for f in filelist if os.path.isfile(f)]

            if len(files_exist) != len(sample_id):
                # create adata file names
                if read_raw_matrix:
                    # create adata file name
                    adata_file_name = os.path.join(saving_adata_dir, 'raw_adata_{}_{}'.format(project_name,
                                                                                              sample_id[ind]))
                else:
                    # create adata file name
                    adata_file_name = os.path.join(saving_adata_dir, 'adata_{}_{}'.format(project_name,
                                                                                          sample_id[ind]))

                print("\nCapture area in Visium slide contains a grid of 4,992 capture spots")

                # 1. Load data
                raw_adata, filtered_adata, config_paths, list_sample_ids = loading_matrices.main(
                    filename=config_filename, read_raw_matrix=read_raw_matrix, spatial_concat=spatial_concat)

                print("-------- Finished: Read out values --------")

                if read_raw_matrix:
                    # save unpre-processed adata object
                    sc.write('{}_unpp.h5'.format(adata_file_name), raw_adata)
                    dataset_processing(adata=raw_adata, results_save_folder=save_folder_ss_analysis,
                                       adata_filename=adata_file_name, filename=config_filename,
                                       single_sample=analysis_sinlge_sample,
                                       preprosessing=pre_processing, dge=dge_analysis, downstream=downstream_analyis,
                                       spatial_velocyto=spatial_velocyto_analysis, filelist=filelist)
                else:
                    # save unpre-processed adata object
                    sc.write('{}_unpp.h5'.format(adata_file_name), filtered_adata)
                    dataset_processing(adata=filtered_adata, results_save_folder=save_folder_ss_analysis,
                                       adata_filename=adata_file_name, filename=config_filename,
                                       single_sample=analysis_sinlge_sample,
                                       preprosessing=pre_processing, dge=dge_analysis, downstream=downstream_analyis,
                                       spatial_velocyto=spatial_velocyto_analysis, filelist=filelist)

            else:
                # adatas exists --> load from storage
                for adata_path_filename in filelist:
                    if read_raw_matrix:
                        raw_adata_tmp = sc.read(adata_path_filename)

                        dataset_processing(adata=raw_adata_tmp, filelist=filelist,
                                           results_save_folder=save_folder_ss_analysis,
                                           adata_filename=adata_path_filename, filename=config_filename,
                                           single_sample=analysis_sinlge_sample,
                                           preprosessing=pre_processing, dge=dge_analysis,
                                           downstream=downstream_analyis, spatial_velocyto=spatial_velocyto_analysis)
                    else:
                        adata_tmp = sc.read(adata_path_filename)
                        dataset_processing(adata=adata_tmp, filelist=filelist,
                                           results_save_folder=save_folder_ss_analysis,
                                           adata_filename=adata_path_filename, filename=config_filename,
                                           single_sample=analysis_sinlge_sample,
                                           preprosessing=pre_processing, dge=dge_analysis,
                                           downstream=downstream_analyis, spatial_velocyto=spatial_velocyto_analysis)
    else:
        # combined analysis on all samples but apply batch correction to correct for technical variability!
        save_folder_analysis = os.path.join(results_save_folder, '_'.join([config_filename]))
        os.makedirs(save_folder_analysis, exist_ok=True)
        if read_raw_matrix:
            # create adata file name
            adata_file_name = os.path.join(saving_adata_dir, 'raw_adata_{}'.format('_'.join([config_filename])))
        else:
            # create adata file name
            adata_file_name = os.path.join(saving_adata_dir, 'adata_{}'.format('_'.join([config_filename])))

        # get files in adata storage path if there are non than file list will be empty []
        filelist = glob.glob(os.path.join(saving_adata_dir, 'adata_{}*.h5'.format('_'.join([config_filename]))))

        # check if pre-processed file exist
        files_exist = [f for f in filelist if '_pp' in f]
        # files_exist = [f for f in filelist if os.path.isfile(f)]

        if len(files_exist) < len([config_filename]):
            print("\nCapture area in Visium slide contains a grid of 4,992 capture spots")
            # 1. Load data
            print("#   --  >Load data and information<  --   #")
            raw_adata, filtered_adata, config_paths, list_sample_ids = loading_matrices.main(
                filename=config_filename, read_raw_matrix=read_raw_matrix, spatial_concat=spatial_concat)
            print("-------- Finished: Read out values --------")

            if read_raw_matrix:
                # save unpre-processed adata object
                sc.write('{}_unpp.h5'.format(adata_file_name), raw_adata)
                dataset_processing(adata=raw_adata, filelist=filelist, results_save_folder=save_folder_analysis,
                                   adata_filename=adata_file_name, filename=config_filename,
                                   single_sample=analysis_sinlge_sample,
                                   preprosessing=pre_processing, dge=dge_analysis, downstream=downstream_analyis,
                                   spatial_velocyto=spatial_velocyto_analysis)
            else:
                # save unpre-processed adata object
                sc.write('{}_unpp.h5'.format(adata_file_name), filtered_adata)
                dataset_processing(adata=filtered_adata, filelist=filelist, results_save_folder=save_folder_analysis,
                                   adata_filename=adata_file_name, filename=config_filename,
                                   single_sample=analysis_sinlge_sample,
                                   preprosessing=pre_processing, dge=dge_analysis, downstream=downstream_analyis,
                                   spatial_velocyto=spatial_velocyto_analysis)
                save_gene_symbol_ensembl_map(filtered_adata, save_folder_analysis)
        else:
            # adata exists --> load from storage
            if pre_processing is True:
                adata_file = [f for f in filelist if '{}_unpp.h5'.format(adata_file_name) in f]
                adata = sc.read(adata_file[0])
            else:
                adata = sc.read(files_exist[0])

            # save_gene_symbol_ensembl_map(adata, save_folder_analysis)

            dataset_processing(adata=adata, filelist=filelist,
                               results_save_folder=save_folder_analysis, adata_filename=adata_file_name,
                               filename=config_filename, single_sample=analysis_sinlge_sample,
                               preprosessing=pre_processing, dge=dge_analysis, downstream=downstream_analyis,
                               spatial_velocyto=spatial_velocyto_analysis)

            # # save adata to csv file in form gene expressions x samples
            # use unpp and QC adata object
            save_count_matrix_to_csv(adata, save_folder_analysis)


if __name__ == '__main__':
    main()
