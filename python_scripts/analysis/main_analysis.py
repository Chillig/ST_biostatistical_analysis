
from python_scripts.analysis.figure_1 import UMAP_bulkRNAseq
from python_scripts.analysis.figure_2 import Fig2A__ST_sections, Fig2B__cytokine_tissuelayer_expression, \
    Fig2C__SC_cytokine_expression, Fig2DF__Bulk_cytokine_expression, Fig2NQ__ST_cytokine_counts
from python_scripts.analysis.figure_3 import Fig3A__ST_UMAP_diseases, Fig3B__ST_UMAP_IL17A_tissuelayers, \
    Fig3CD__ST_Volcanoplot_cytokine
from python_scripts.analysis.figure_4 import Fig4A__SC_UMAP_clusters, Fig4B__SC_UMAP_IL17A, Fig4C__SC_Volcano_plot
from python_scripts.analysis.figure_5 import Fig5AC__ST_pseudobulk_aggregation_Correlation, \
    Fig4FH__Weighted_Correlation
# Creates Figure 5G-I
from python_scripts.spatial_correlation import main as csdcc

from python_scripts.analysis.supplfigure_1 import SuppFig1C__cytokine_counts_skinlayers, \
    SuppFig1DE__Lesion_NonLesion_Mean_expression_cytokines_others
from python_scripts.analysis.supplfigure_3 import SuppFig3ABC__ST_UMAP
from python_scripts.analysis.supplfigure_4 import SupplFig4AF__spot_celltypes
from python_scripts.analysis.supplfigure_5 import SuppFig5AC__ST_Volcano_Boxplot
from python_scripts.analysis.supplfigure_6 import SuppFig6A__SC_UMAP_IFNG, SuppFig6B__SC_Volcano_IFNG
from python_scripts.analysis.supplfigure_7 import SuppFig7H__cytokines_vs_responders

from python_scripts.analysis.suppltable import SuppTab2__overview_counts

from python_scripts.utils import gene_lists

import os
from datetime import date
import scanpy as sc
import numpy as np
import pandas as pd
from operator import itemgetter

import matplotlib.pyplot as plt


def get_deg_files(input_dir, cytokine):
    design_function_spatial = 'cdr_project_patient_annotation_cyto'
    design_function_singlecell = 'cdr_annotation_cyto'
    dataset = "Whole_T_cell_matrix"

    # load df's using pandas
    # Spatial DEGs paths
    comparisons = "_".join([cytokine, 'vs_Others'])
    st_path = os.path.join(
        input_dir, "dge_analysis",
        "".join(['2022-04-08_spatial_', design_function_spatial, os.path.sep, 'spatial',
                 os.path.sep, dataset, "__", design_function_spatial, os.path.sep, cytokine, os.path.sep,
                 comparisons]))
    st_file = os.path.join(
        input_dir, "dge_analysis", st_path,
        "".join([dataset, "_", comparisons, "_glmGamPoi_DGE_all_genes.csv"]))

    # Single cell DEGs paths
    sc_path = os.path.join(
        input_dir, "dge_analysis", '2021-02-01_single_cell__', design_function_singlecell, 'single_cell',
        "".join([dataset, "__", design_function_singlecell]), cytokine, comparisons)
    sc_file = os.path.join(
        input_dir, "dge_analysis", sc_path,
        "".join(["single_cell_", comparisons, "_glmGamPoi_DGE_all_genes.csv"]))

    return st_path, st_file, sc_path, sc_file


def load_spatial_data(adata_folder):
    date_st_unpp = '2022-04-08'
    date_st_pp = '2022-04-08'

    # load unpreprocessed data
    unpp_adata = sc.read(os.path.join(
        adata_folder, date_st_unpp, "Spatial Transcriptomics_unpp_cleaned_LPADPso.h5"))

    # load preprocessed data
    pp_adata = sc.read(os.path.join(adata_folder, date_st_pp, 'st_QC_normed_BC_project_PsoADLP.h5'))

    return unpp_adata, pp_adata


def load_singlecell_data(adata_folder):
    unpp_adata = sc.read(os.path.join(adata_folder, 'single_cell', '2020-11-30', 'sc_adata_unpp.h5'))
    pp_adata = sc.read(os.path.join(adata_folder, '2020-12-04_SC_Data_QC_clustered.h5'))

    return unpp_adata, pp_adata


def load_bulk_data(data_folder):
    # --> Bulk RNAseq data
    input_rnaseq = os.path.join(data_folder, "bulk_RNAseq")
    # - Read bulk-RNAseq count matrix
    data = pd.read_csv(os.path.join(input_rnaseq, "bulkRNA_countMat.txt"), sep='\t')
    # - Read in metaData
    metadata = pd.read_excel(os.path.join(input_rnaseq, "bulkRNA_metaData.xlsx"))

    return data, metadata


def main(save_folder: str, adata_folder: str, input_folder: str):
    # Data
    # --> Spatial Transcriptomics data
    unpp_st_adata, pp_st_adata = load_spatial_data(adata_folder)
    # --> Single cell RNAseq data
    sc_adata_folder = '/Users/christina.hillig/Documents/Projects/annData_objects'
    unpp_sc_adata, pp_sc_adata = load_singlecell_data(adata_folder=sc_adata_folder)
    # --> Bulk RNAseq data
    bulk_data, meta_data = load_bulk_data(data_folder=input_folder)

    # sample: technical sample
    # project: from sequencing batch the samples originate
    # capture area: a capture area on the Visium object slide
    # object slide: Visium object slide
    # specimen: a tissue slice
    # SAMPLE: clinical patient IDs
    # patients: patient IDs 1-40

    """ Input parameter """
    st_path_il17a, st_file_il17a, sc_path_il17a, sc_file_il17a = get_deg_files(input_dir=input_folder, cytokine='IL17A')
    st_path_ifng, st_file_ifng, sc_path_ifng, sc_file_ifng = get_deg_files(input_dir=input_folder, cytokine='IFNG')
    input_st_dgeanalysis = os.path.join(
        input_folder, 'dge_analysis', '2022-04-08_spatial__cdr_project_patient_annotation_cyto')

    # Spatial Correlation
    radius = list(np.arange(0, 10, 1))
    corr_method = 'spearman'
    t_cell_cytocines, cyto_resps_list, cytokine_responders = gene_lists.get_publication_cyto_resps()

    """ ----  Main Figures  ---- """
    """ Figure 1 """
    save_folder_fig1 = os.path.join(save_folder, 'Fig1', str(date.today()))
    os.makedirs(save_folder_fig1, exist_ok=True)

    UMAP_bulkRNAseq.main(save_folder=save_folder_fig1, bulk_rnaseq=bulk_data, metadata=meta_data)

    """ Figure 2 """
    Fig2A__ST_sections.main(save_folder=save_folder, adata=unpp_st_adata)
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
    save_folder_fig3 = os.path.join(save_folder, 'Fig3', str(date.today()))
    os.makedirs(save_folder_fig3, exist_ok=True)

    Fig3A__ST_UMAP_diseases.main(save_folder=save_folder_fig3, spatial_adata=pp_st_adata)
    Fig3B__ST_UMAP_IL17A_tissuelayers.main(save_folder=save_folder_fig3, spatial_adata=pp_st_adata)
    # Figure 3CD: Plots Volcano plot (C) and Violinplots (D)
    # df_keys contains column names of DEGs dataframe.
    # First string should be the name of the Log2FC column,
    # Second of either p-value or p-adjusted value,
    # Third should be the name of the HNGO gene names
    Fig3CD__ST_Volcanoplot_cytokine.main(
        adata=pp_st_adata, save_folder=save_folder_fig3, df_keys=['log2fc', 'padj', 'gene_symbol'],
        log=False, dge_results_folder=input_st_dgeanalysis)

    """ Figure 4 """
    save_folder_fig4 = os.path.join(save_folder, 'Fig4', str(date.today()))
    os.makedirs(save_folder_fig4, exist_ok=True)

    Fig4A__SC_UMAP_clusters.main(save_folder=save_folder, pp_adata=pp_sc_adata, cluster_algorithm='leiden')
    Fig4B__SC_UMAP_IL17A.main(save_folder=save_folder, adata=pp_sc_adata)
    # Figure 4C
    # df_keys contains column names of DEGs dataframe.
    # First string should be the name of the Log2FC column,
    # Second of either p-value or p-adjusted value,
    # Third should be the name of the HNGO gene names
    Fig4C__SC_Volcano_plot.main(
        dataset_type='SC', save_folder=save_folder_fig4, df_keys=['log2fc', 'padj', 'gene_symbol'],
        log=False, dge_results_folder=sc_path_il17a)

    """ Figure 5 """
    save_folder_fig5 = '/Volumes/CH__data/ST_immune_publication/Revision/Fig5'
    save_folder_fig5ac = os.path.join(save_folder_fig5, 'Bulk_Weighted_Spearman', 'Unweighted_fit', str(date.today()))
    os.makedirs(save_folder_fig5ac, exist_ok=True)
    # Run Weighted correlation by cytokine+ spots in epidermis for pseudo-bulk approach
    fig5__dict_weighted_transcripts_corr, df_bulk = Fig5AC__ST_pseudobulk_aggregation_Correlation.main(
        save_folder=save_folder_fig5ac, adata=unpp_st_adata, corr_method=corr_method)
    # Read out correlation values for each cytokine
    data = list(itemgetter(*[0, 4, 8])(fig5__dict_weighted_transcripts_corr['spearman']))
    df_fig5ac = pd.DataFrame(data, columns=['Cytokine', 'Responders', 'Spearman Correlation', 'P-value'])

    # Create figure 5D-H - if get_plots is set to True, otherwise it will only create the plots for 5E-H
    save_folder_fig5eg = os.path.join(save_folder, 'Weighted_Spearman_unppadata', 'Unweighted_fit', str(date.today()))
    os.makedirs(save_folder_fig5eg, exist_ok=True)
    adata, counts_dict, df_stats_responders_in_vs_outlesion_sdc, \
    df_stats_cytokines_responders_in_sdc, df_radius_vs_spearman = csdcc.main(
        save_folder=save_folder_fig5eg, adata=unpp_st_adata, corr_method=corr_method, get_plots=False,
        find_associated_genes=False, radius=radius,
        cond_genes=t_cell_cytocines, genes_resps=cytokine_responders)
    max_corr_value = df_radius_vs_spearman[
        df_radius_vs_spearman.columns[df_radius_vs_spearman.columns.str.contains('correlation')]].idxmax(axis=0)

    # Read out counts and correlation values
    list_fig5eh_counts = list(itemgetter(*max_corr_value)(counts_dict))
    df_fig5eh_counts = pd.concat(
        [pd.DataFrame(list_fig5eh_counts[0]), pd.DataFrame(list_fig5eh_counts[1]), pd.DataFrame(list_fig5eh_counts[2])],
        axis=0, keys=[str(max_corr_value.values[0]), str(max_corr_value.values[1]), str(max_corr_value.values[2])])
    df_fig5eh_correlation = df_radius_vs_spearman.iloc[max_corr_value]

    writer = pd.ExcelWriter(os.path.join(save_folder_fig5eg, 'Figure_5.xlsx'), engine='xlsxwriter')
    df_fig5ac.to_excel(writer, sheet_name='AC_correlation')
    df_bulk.to_excel(writer, sheet_name='AC_counts')
    df_radius_vs_spearman.to_excel(writer, sheet_name='E__radius_vs_correlation')
    df_fig5eh_correlation.to_excel(writer, sheet_name='FH_correlation')
    df_fig5eh_counts.to_excel(writer, sheet_name='FH_counts')
    writer.save()
    writer.close()

    """ --- Workflow Figure 1EF --- """
    # Have to run Correlation analysis before
    # H&E image for workflow Fig 1
    adata.obs['IL17A_Responders_label'] = adata.obs['IL17A_Responders_clusters'].copy()
    adata.obs['IL17A_Responders_label'] = adata.obs['IL17A_Responders_label'].astype(str)
    adata.obs['IL17A_Responders_label'][adata.obs['IL17A_clusters'] == 'IL17A'] = 'IL17A'
    adata.obs['IL17A_Responders_label'] = adata.obs['IL17A_Responders_label'].astype('category')
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.spatial(adata[adata.obs['sample'] == 'P15509_1004'], color='IL17A_Responders_label', size=1,
                  library_id='P15509_1004', ax=ax, palette=['#377eb8', '#B8860B'],
                  add_outline=True, alpha=1, outline_width=(8, 0.05), show=False)
    plt.tight_layout()
    fig.savefig(os.path.join(save_folder, 'Workflow_Fig1E.pdf'))
    plt.close(fig=fig)

    adata.obs['IL17A_in_sdcc_r1'] = adata.obs['IL17A_in_sdcc_r1'].astype('category')
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.spatial(adata[(adata.obs['sample'] == 'P15509_1004') & (adata.obs['IL17A_in_sdcc_r1'] != 0)],
                  color='IL17A_in_sdcc_r1', size=1, library_id='P15509_1004', ax=ax, palette=['#377eb8', '#B8860B'],
                  add_outline=True, alpha=1, outline_width=(10, 0.05), show=False)
    plt.tight_layout()
    fig.savefig(os.path.join(save_folder, 'Workflow_Fig1F.pdf'))
    plt.close(fig=fig)

    # ------------------------------------

    """ ----  Supplemental Figures  ---- """
    """ Figure S1: NcISD are characterized by low cytokine UMI counts in skin """
    save_folder_sfig1 = os.path.join(save_folder, 'SFig1')
    os.makedirs(save_folder_sfig1, exist_ok=True)
    SuppFig1C__cytokine_counts_skinlayers.main(save_folder=save_folder, adata=unpp_st_adata)
    df_stats_mean_cytokines_others_l_vs_nl = SuppFig1DE__Lesion_NonLesion_Mean_expression_cytokines_others.main(
        save_folder=save_folder_sfig1, pp_st_adata=pp_st_adata)

    """ Figure S3: Cytokine transcript-positive spots are located in the epidermis and enriched in lesional skin """
    save_folder_sfig3 = os.path.join(save_folder, 'SFig3')
    os.makedirs(save_folder_sfig3, exist_ok=True)
    SuppFig3ABC__ST_UMAP.main(
        save_folder=save_folder_sfig3, spatial_adata=pp_st_adata, spatial_cluster_label='tissue_layer')

    """ Figure S4: Tangram analysis of cytokine-transcript positive spatial spots reveals heterogeneous
    cellular composition """
    save_folder_sfig4 = os.path.join(save_folder, 'SFig4')
    os.makedirs(save_folder_sfig4, exist_ok=True)
    print("Run jupyter-notebooks before you call this function")
    SupplFig4AF__spot_celltypes.main(save_folder=save_folder_sfig3)

    """ Figure S5: Cytokine transcript-positive leukocyte spots are characterized by specific gene signatures """
    save_folder_sfig5 = os.path.join(save_folder, 'SFig5')
    os.makedirs(save_folder_sfig5, exist_ok=True)
    SuppFig5AC__ST_Volcano_Boxplot.main(
        adata=pp_st_adata, save_folder=save_folder_sfig4, df_keys=['log2fc', 'padj', 'gene_symbol'],
        log=False, dge_results_folder=st_path_ifng)

    """ Figure S6: Single-cell analysis reveals specific gene signatures for IL17A and IFNG expressing cells """
    save_folder_sfig6 = os.path.join(save_folder, 'SFig6')
    os.makedirs(save_folder_sfig6, exist_ok=True)
    SuppFig6A__SC_UMAP_IFNG.main(save_folder=save_folder, adata=pp_sc_adata)
    SuppFig6B__SC_Volcano_IFNG.main(
        dataset_type='SC', save_folder=save_folder_sfig6, df_keys=['log2fc', 'padj', 'gene_symbol'],
        log=False, dge_results_folder=sc_path_ifng)

    """ Figure S7H: Identification of cytokine responder genes for spatial correlation """
    save_folder_sfig7 = os.path.join(save_folder, 'SFig7')
    os.makedirs(save_folder_sfig7, exist_ok=True)
    # For Figures E-G see r_scripts
    df_stats_lesion_mean_persample_perspecimen_cytokines_responders = SuppFig7H__cytokines_vs_responders.main(
        save_folder=save_folder_sfig7, pp_st_adata=pp_st_adata)

    """ Figure S8: Spatial transcriptomics enables to identify new cytokine-related responder genes and thereby
    potential drug targets """
    save_folder_sfig8 = os.path.join(save_folder, 'Figure_S8', str(date.today()))
    os.makedirs(save_folder_sfig8, exist_ok=True)
    # 1. Run with preprocessed (QC, normed, bc) data density clustering
    _, _, _ = csdcc.main(save_folder=save_folder_sfig8, adata=pp_st_adata, corr_method=corr_method, get_plots=False,
                         find_associated_genes=True, radius=list(np.arange(0, 5, 1)), cond_genes=t_cell_cytocines,
                         genes_resps=cytokine_responders)
    # 2. Identify DEGs (python) and create Volcano plots (R script)
    # 3. Identify cytokine associated genes using DEGs from optimal radii 4, 3, and 0 for IFNG, IL13 and IL17A

    """ Supplemental Tables """
    # Supplemental Table 2
    SuppTab2__overview_counts.create_supptable_1(unpp_st_adata=unpp_st_adata, save_folder=save_folder)

    # Supplemental Table 3
    writer = pd.ExcelWriter(os.path.join(save_folder, 'Supplemental_Table_3.xlsx'), engine='xlsxwriter')
    df_stats_responders_in_vs_outlesion_sdc.to_excel(writer, sheet_name='Abstract')
    df_stats_mean_cytokines_others_l_vs_nl.to_excel(writer, sheet_name='FigS1A&B')
    df_stats_lesion_mean_persample_perspecimen_cytokines_responders.to_excel(writer, sheet_name='FigS7H')
    writer.save()
    writer.close()


if __name__ == '__main__':
    today = date.today()
    project_folder = os.path.join("..", "..")
    output_path = os.path.join(project_folder, "output", "reviewers")
    os.makedirs(output_path, exist_ok=True)
    adata_path = os.path.join(project_folder, "adata_storage")
    input_path = os.path.join(project_folder, "input")
    main(save_folder=output_path, adata_folder=adata_path, input_folder=input_path)
