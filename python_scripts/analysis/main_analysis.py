
from python_scripts.analysis.figure_1 import UMAP_bulkRNAseq
from python_scripts.analysis.figure_2 import Fig2A__ST_sections, Fig2B__cytokine_tissuelayer_expression, \
    Fig2C__SC_cytokine_expression, Fig2DF__Bulk_cytokine_expression, Fig2NQ__ST_cytokine_counts
from python_scripts.analysis.figure_3 import Fig3A__ST_UMAP_IL17A_tissuelayers, Fig3B__ST_Volcanoplot_cytokine
from python_scripts.analysis.figure_4 import Fig4A__SC_UMAP_clusters, Fig4B__SC_UMAP_IL17A, Fig4C__SC_Volcano_plot
from python_scripts.analysis.figure_5 import Fig5AC__ST_pseudobulk_aggregation_Correlation, \
    Fig4FH__Weighted_Correlation
# Creates Figure 5G-I
from python_scripts.spatial_correlation import main as csdcc

from python_scripts.analysis.supplfigure_1 import SuppFig1A__cytokine_counts_skinlayers
from python_scripts.analysis.supplfigure_2 import SuppFig2ABC__ST_UMAP
from python_scripts.analysis.supplfigure_3 import SupplFig3AF__spot_celltypes
from python_scripts.analysis.supplfigure_4 import SuppFig4AC__ST_Volcano_Boxplot
from python_scripts.analysis.supplfigure_5 import SuppFig5A__SC_UMAP_IFNG, SuppFig5B__SC_Volcano_IFNG
from python_scripts.analysis.supplfigure_7 import SuppFig7AC__ST_Common_genes_cytokinepositive_condition
from python_scripts.analysis.supplfigure_8 import SuppFig8A__ST_pseudobulk_Correlation_permutedresponders, \
    SuppFig8AC__refined_responds

from python_scripts.analysis.suppltable import SuppTab1__overview_counts

from python_scripts.utils import gene_lists

import os
from datetime import date
import scanpy as sc
import numpy as np
import pandas as pd


def get_deg_files(input_dir, cytokine):
    design_function = 'cdr_project_patient_annotation_cyto'
    dataset = "Whole_T_cell_matrix"

    # load df's using pandas
    comparisons = "_".join([cytokine, 'vs_Others'])
    st_path = os.path.join(
        input_dir, "dge_analysis",
        "".join(['2022-04-08_spatial__cdr_project_patient_annotation_cyto', os.path.sep, 'spatial',
                 os.path.sep, dataset, "__", design_function, os.path.sep, cytokine, os.path.sep,
                 comparisons]))
    st_file = os.path.join(
        input_dir, "dge_analysis", st_path,
        "".join(["Whole_T_cell_matrix_", comparisons, "_glmGamPoi_DGE_all_genes.csv"]))

    sc_path = os.path.join(
        input_dir, "dge_analysis", '2021-02-01_single_cell__cdr_annotation_cyto', 'single_cell',
        "".join([dataset, "__", "cdr_annotation_cyto"]), cytokine, comparisons)
    sc_file = os.path.join(
        input_dir, "dge_analysis", sc_path,
        "".join(["single_cell_", comparisons, "_glmGamPoi_DGE_all_genes.csv"]))

    return st_path, st_file, sc_path, sc_file


def main(save_folder: str, adata_folder: str, input_folder: str):
    date_st_unpp = '2022-04-08'  # "2020-10-06" -> "st_adata_P15509_P16357_wo_4_7_unpp.h5"
    date_st_pp = '2022-04-08'  # "2020-12-04" -> 2020-12-04_Visium_Data_QC_BC_clustered.h5
    unpp_st_adata = sc.read(os.path.join(adata_folder, date_st_unpp, "Spatial Transcriptomics_unpp_cleaned.h5"))
    unpp_st_adata.obs.loc[(unpp_st_adata.obs['basal EPIDERMIS'] == 1) & (unpp_st_adata.obs['DERdepth1'] == 1),
                          'basal EPIDERMIS'] = [0, 0, 1]
    unpp_st_adata.obs.loc[(unpp_st_adata.obs['basal EPIDERMIS'] == 1) & (unpp_st_adata.obs['DERdepth1'] == 1),
                          'DERdepth1'] = 0
    unpp_st_adata = unpp_st_adata[unpp_st_adata.obs['DISEASE'] != 'PRP'].copy()
    pp_st_adata = sc.read(os.path.join(adata_folder, date_st_pp, 'st_QC_normed_BC_project_PsoADLP.h5'))
    pp_st_adata.obs.loc[(pp_st_adata.obs['basal EPIDERMIS'] == 1) & (pp_st_adata.obs['DERdepth1'] == 1),
                        'basal EPIDERMIS'] = [0, 0, 1]
    pp_st_adata.obs.loc[(pp_st_adata.obs['basal EPIDERMIS'] == 1) & (pp_st_adata.obs['DERdepth1'] == 1),
                        'DERdepth1'] = 0

    # sample: technical sample
    # project: from sequencing batch the samples originate
    # capture area: a capture area on the Visium object slide
    # object slide: Visium object slide
    # specimen: a tissue slice
    # SAMPLE: clinical patient IDs
    # patients: patient IDs 1-40

    # TODO check if number of samples influences the results -> umap do they cluster together?
    sc_adata_folder = '/Users/christina.hillig/Documents/Projects/annData_objects'
    unpp_sc_adata = sc.read(os.path.join(sc_adata_folder, '2020-11-30', 'sc_adata_unpp.h5'))
    pp_sc_adata = sc.read(os.path.join(sc_adata_folder, '2020-12-04_SC_Data_QC_clustered.h5'))
    # --> Bulk RNAseq data
    input_rnaseq = os.path.join(input_folder, "bulk_RNAseq")
    # - Read bulk-RNAseq count matrix
    bulk_data = pd.read_csv(os.path.join(input_rnaseq, "bulkRNA_countMat.txt"), sep='\t')
    # - Read in metaData
    meta_data = pd.read_excel(os.path.join(input_rnaseq, "bulkRNA_metaData.xlsx"))

    """ Input parameter """
    st_path_il17a, st_file_il17a, sc_path_il17a, sc_file_il17a = get_deg_files(input_dir=input_folder, cytokine='IL17A')
    st_path_ifng, st_file_ifng, sc_path_ifng, sc_file_ifng = get_deg_files(input_dir=input_folder, cytokine='IFNG')
    # Spatial Correlation
    radius = list(np.arange(0, 10, 1))
    corr_method = 'spearman'
    sfig3f_input_path = os.path.join(
        save_folder, "Whole_T_cell_matrix__cdr_project_patient_annotation_cyto")

    """ ----  Main Figures  ---- """
    """ Figure 1 """
    UMAP_bulkRNAseq.main(save_folder=save_folder, bulk_rnaseq=bulk_data, metadata=meta_data)

    """ Figure 2 """
    Fig2A__ST_sections.main(save_folder=save_folder, adata=unpp_st_adata)
    # TODO check if you can apply wilcoxon test ..
    Fig2B__cytokine_tissuelayer_expression.main(save_folder=save_folder, adata=unpp_st_adata)
    Fig2C__SC_cytokine_expression.main(save_folder=save_folder, adata=unpp_sc_adata)
    Fig2DF__Bulk_cytokine_expression.main(save_folder=save_folder, bulk_rnaseq=bulk_data, metadata=meta_data)
    # Fig2N-P: Cytokine counts of IFNG, IL17A, IL13 resolved by disease
    Fig2NQ__ST_cytokine_counts.get_per_diagnosis_cytokine_counts(unpp_st_adata, save_folder, biopsy_type='LESIONAL')
    Fig2NQ__ST_cytokine_counts.get_per_diagnosis_cytokine_counts(unpp_st_adata, save_folder, biopsy_type='NON LESIONAL')
    # Fig2Q: Total number of cytokine counts of IFNG, IL17A, IL13 resolved by disease
    Fig2NQ__ST_cytokine_counts.get_cytokinecounts_per_diagnosis(unpp_st_adata, save_folder, biopsy_type='LESIONAL')
    Fig2NQ__ST_cytokine_counts.get_cytokinecounts_per_diagnosis(unpp_st_adata, save_folder, biopsy_type='NON LESIONAL')

    """ Figure 3 """
    Fig3A__ST_UMAP_IL17A_tissuelayers.main(save_folder=save_folder, spatial_adata=pp_st_adata)
    input_st_dgeanalysis = os.path.join(
        input_folder, 'dge_analysis', '2022-04-08_spatial__cdr_project_patient_annotation_cyto')
    Fig3B__ST_Volcanoplot_cytokine.main(
        adata=pp_st_adata, save_folder=save_folder, df_keys=['log2fc', 'pval', 'gene_symbol'],
        log=False, dge_results_folder=input_st_dgeanalysis)

    """ Figure 4 """
    Fig4A__SC_UMAP_clusters.main(save_folder=save_folder, pp_adata=pp_sc_adata, cluster_algorithm='leiden')
    Fig4B__SC_UMAP_IL17A.main(save_folder=save_folder, adata=pp_sc_adata)
    Fig4C__SC_Volcano_plot.main(dataset_type='SC', save_folder=save_folder, df_keys=['log2fc', 'pval', 'gene_symbol'],
                                log=False, dge_results_folder=sc_path_il17a)

    """ Figure 5 """
    t_cell_cytocines, cyto_resps_list, cytokine_responders = gene_lists.get_publication_cyto_resps()
    save_folder = '/Volumes/CH__data/ST_immune_publication/Revision/Fig5'
    save_folder_fig5ac = os.path.join(save_folder, 'Bulk_Weighted_Spearman', 'Unweighted_fit', str(date.today()))
    os.makedirs(save_folder_fig5ac, exist_ok=True)
    # Run Weighted correlation by cytokine+ spots in epidermis for pseudo-bulk approach
    fig5__dict_weighted_transcripts_corr, df_bulk = Fig5AC__ST_pseudobulk_aggregation_Correlation.main(
            save_folder=save_folder_fig5ac, adata=unpp_st_adata, corr_method=corr_method)

    # Create figure 5E-G
    save_folder_fig5eg = os.path.join(save_folder, 'Weighted_Spearman_unppadata', 'Unweighted_fit', str(date.today()))
    os.makedirs(save_folder_fig5eg, exist_ok=True)
    counts_dict = csdcc.main(
        save_folder=save_folder_fig5eg, adata=unpp_st_adata, corr_method=corr_method, get_plots=False,
        find_responders=False, radius=radius,
        cond_genes=t_cell_cytocines, genes_resps=cytokine_responders)

    """ ----  Supplemental Figures  ---- """
    SuppFig1A__cytokine_counts_skinlayers.main(save_folder=save_folder, adata=unpp_st_adata)

    """ SFig 2 """
    SuppFig2ABC__ST_UMAP.main(save_folder=save_folder, spatial_adata=pp_st_adata, spatial_cluster_label='tissue_layer')

    """ SFig 3 """  # TODO add SFig3 Tangram jupyternotebook..
    SupplFig3AF__spot_celltypes.main(save_folder=save_folder)

    """ SFig 4 """
    SuppFig4AC__ST_Volcano_Boxplot.main(
        adata=pp_st_adata, save_folder=save_folder, df_keys=['log2fc', 'pval', 'gene_symbol'],
        log=False, dge_results_folder=st_path_ifng)

    """ SFig 5 """
    SuppFig5A__SC_UMAP_IFNG.main(save_folder=save_folder, adata=pp_sc_adata)
    SuppFig5B__SC_Volcano_IFNG.main(
        dataset_type='SC', save_folder=save_folder, df_keys=['log2fc', 'pval', 'gene_symbol'],
        log=False, dge_results_folder=sc_path_ifng)

    """ SFig 6"""

    """ SFig 7 """
    SuppFig7AC__ST_Common_genes_cytokinepositive_condition.main(input_path=st_path_ifng, save_folder=save_folder)

    """ SFig 8"""
    save_folder_sfig8 = os.path.join(save_folder, 'Figure_10', str(date.today()))
    os.makedirs(save_folder_sfig8, exist_ok=True)
    # 1. Run with preprocessed (QC, normed, bc) data density clustering
    csdcc.main(save_folder=save_folder_sfig8, adata=pp_st_adata, corr_method=corr_method, get_plots=False,
               find_responders=True, radius=list(np.arange(0, 6, 1)), cond_genes=t_cell_cytocines,
               genes_resps=cytokine_responders)
    # 2. Identify DEGs (python) and create Volcano plots (R script)
    # 3. Identify new responder genes using DEGs from radii 1-5 and rerun density clustering
    SuppFig8AC__refined_responds.main()

    """ Supplemental Table """
    # Supplemental Table 1
    SuppTab1__overview_counts.create_supptable_1(unpp_st_adata=unpp_st_adata, save_folder=save_folder)


if __name__ == '__main__':
    today = date.today()
    project_folder = os.path.join("..", "..")
    output_path = os.path.join(project_folder, "output", "reviewers")
    os.makedirs(output_path, exist_ok=True)
    adata_path = os.path.join(project_folder, "adata_storage")
    input_path = os.path.join(project_folder, "input")
    main(save_folder=output_path, adata_folder=adata_path, input_folder=input_path)
