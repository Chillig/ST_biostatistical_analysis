
from python_scripts.analysis.figure_1 import UMAP_bulkRNAseq
from python_scripts.analysis.figure_2 import Fig2A__ST_sections, Fig2B__cytokine_tissuelayer_expression, \
    Fig2C__SC_cytokine_expression, Fig2DF__Bulk_cytokine_expression
from python_scripts.analysis.figure_3 import Fig3A__ST_UMAP_IL17A_tissuelayers, Fig3B__ST_Volcanoplot_cytokine, \
    Fig3D__SC_UMAP_IL17A, Fig3E__STSC_Correlation_IL17A
from python_scripts.analysis.figure_4 import Fig4A__ST_pseudobulk_aggregation_Correlation, \
    Fig4B__Image_CytokineResponders
# Creates Figure 4C
from python_scripts.spatial_correlation import main as csdcc

from python_scripts.analysis.supplfigure_1 import SuppFig1A__cytokine_counts_skinlayers
from python_scripts.analysis.supplfigure_2 import SuppFig2ABC__ST_UMAP
from python_scripts.analysis.supplfigure_3 import SuppFig3AC__ST_UMAP_IL13_IFNG_skinlayers, \
    SuppFig3E__STSC_Correlation_IFNG, SuppFig3F__ST_Common_genes_cytokinepositive_condition
from python_scripts.analysis.supplfigure_4 import SuppFig4A__SC_UMAP_clusters, SuppFig4B__SC_UMAP_IFNG, \
    SuppFig4CD__SC_Volcano_cytokines
from python_scripts.analysis.supplfigure_7 import SuppFig7A__ST_pseudobulk_Correlation_permutedresponders


import os
from datetime import date
import scanpy as sc
import pandas as pd


def get_deg_files(input_dir, cytokine):
    design_function = 'cdr_patient_annotation_cyto'
    dataset = "Whole_T_cell_matrix"

    # load df's using pandas
    comparisons = "_".join([cytokine, 'vs_Others'])
    st_path = os.path.join(
        input_dir, "dge_analysis",
        "".join(['2021-02-01_spatial__cdr_patient_annotation_cyto', os.path.sep, 'spatial',
                 os.path.sep, dataset, "__", design_function, os.path.sep, cytokine, os.path.sep,
                 comparisons, os.path.sep, "spatial_", comparisons, "_glmGamPoi_DGE_all_genes.csv"]))

    sc_path = os.path.join(
        input_dir, "dge_analysis",
        "".join(['2021-02-01_single_cell__cdr_annotation_cyto', os.path.sep, 'single_cell', os.path.sep,
                 dataset, "__", "cdr_annotation_cyto", os.path.sep, cytokine,
                 os.path.sep, comparisons, os.path.sep, "single_cell_", comparisons, "_glmGamPoi_DGE_all_genes.csv"]))

    return st_path, sc_path


def main(save_folder, adata_folder, input_folder):
    unpp_st_adata = sc.read(os.path.join(adata_folder, "2020-10-06", "st_adata_P15509_P16357_wo_4_7_unpp.h5"))
    pp_st_adata = sc.read(os.path.join(adata_folder, "2020-12-04_Visium_Data_QC_BC_clustered.h5"))
    unpp_sc_adata = sc.read(os.path.join(adata_folder, '2020-11-30', 'sc_adata_unpp.h5'))
    pp_sc_adata = sc.read(os.path.join(adata_folder, '2020-12-04_SC_Data_QC_clustered.h5'))
    # --> Bulk RNAseq data
    input_rnaseq = os.path.join(input_folder, "bulk_RNAseq")
    # - Read bulk-RNAseq count matrix
    bulk_data = pd.read_csv(os.path.join(input_rnaseq, "bulkRNA_countMat.txt"), sep='\t')
    # - Read in metaData
    meta_data = pd.read_excel(os.path.join(input_rnaseq, "bulkRNA_metaData.xlsx"))

    """ ----  Main Figures  ---- """
    # figure 1 plot
    UMAP_bulkRNAseq.main(save_folder=save_folder, bulk_rnaseq=bulk_data, metadata=meta_data)

    # Figure 2
    Fig2A__ST_sections.main(save_folder=save_folder, adata=unpp_st_adata)
    Fig2B__cytokine_tissuelayer_expression.main(save_folder=save_folder, adata=unpp_st_adata)
    Fig2C__SC_cytokine_expression.main(save_folder=save_folder, adata=unpp_sc_adata)
    Fig2DF__Bulk_cytokine_expression.main(save_folder=save_folder, bulk_rnaseq=bulk_data, metadata=meta_data)

    # Figure 3
    Fig3A__ST_UMAP_IL17A_tissuelayers.main(save_folder=save_folder, spatial_adata=pp_st_adata)
    input_st_dgeanalysis = os.path.join(input_folder, 'dge_analysis', '2021-02-01_spatial__cdr_patient_annotation_cyto')
    Fig3B__ST_Volcanoplot_cytokine.main(
        dataset_type='ST', save_folder=save_folder, df_keys=['log2fc', 'pval', 'gene_symbol'],
        log=False, dge_results_folder=input_st_dgeanalysis)
    Fig3D__SC_UMAP_IL17A.main(save_folder=save_folder, adata=pp_sc_adata)
    st_path, sc_path = get_deg_files(input_dir=input_folder, cytokine='IL17A')
    Fig3E__STSC_Correlation_IL17A.main(path_st_data=st_path, path_sc_data=sc_path, save_folder=save_folder,
                                       log2fc_cut=1., pval_cut=0.05, zoom=True)

    # Figure 4
    Fig4A__ST_pseudobulk_aggregation_Correlation.main(save_folder=save_folder, adata=unpp_st_adata)
    Fig4B__Image_CytokineResponders.main(save_folder=save_folder, adata=unpp_st_adata)
    csdcc.main(save_folder=save_folder, adata=unpp_st_adata)

    """ ----  Supplemental Figures  ---- """
    SuppFig1A__cytokine_counts_skinlayers.main(save_folder=save_folder, adata=unpp_st_adata)
    SuppFig2ABC__ST_UMAP.main(save_folder=save_folder, spatial_adata=pp_st_adata)
    SuppFig3AC__ST_UMAP_IL13_IFNG_skinlayers.main(save_folder=save_folder, spatial_adata=pp_st_adata)
    st_path, sc_path = get_deg_files(input_dir=input_folder, cytokine='IFNG')
    SuppFig3E__STSC_Correlation_IFNG.main(path_st_data=st_path, path_sc_data=sc_path, save_folder=save_folder,
                                          log2fc_cut=1., pval_cut=0.05, zoom=True)
    SuppFig3F__ST_Common_genes_cytokinepositive_condition.main(save_folder=save_folder)
    SuppFig4A__SC_UMAP_clusters.main(save_folder=save_folder, pp_adata=pp_sc_adata, cluster_algorithm='leiden')
    SuppFig4B__SC_UMAP_IFNG.main(save_folder=save_folder, adata=pp_sc_adata)
    input_sc_dgeanalysis = os.path.join(input_folder, 'dge_analysis', '2021-02-01_single_cell__cdr_annotation_cyto')
    SuppFig4CD__SC_Volcano_cytokines.main(
        dataset_type='SC', save_folder=save_folder, df_keys=['log2fc', 'pval', 'gene_symbol'],
        log=False, dge_results_folder=input_sc_dgeanalysis)
    SuppFig7A__ST_pseudobulk_Correlation_permutedresponders.main(save_folder=save_folder, adata=unpp_st_adata)

    # TODO add Figure 4C and Suppl Fig 7B


if __name__ == '__main__':
    today = date.today()
    project_folder = os.path.join("..", "..")
    output_path = os.path.join(project_folder, "output")
    adata_path = os.path.join(project_folder, "adata_storage")
    input_path = os.path.join(project_folder, "input")
    main(save_folder=output_path, adata_folder=adata_path, input_folder=input_path)
