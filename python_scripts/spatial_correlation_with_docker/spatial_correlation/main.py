import os
import numpy as np
import pandas as pd
from datetime import date
import scanpy as sc

# import python scripts
import gene_lists
import corr_statistics
import evaluations
import plot_spatial_correlation
# import classes to run clustering algorithm
import modul_density_clustering


def main():
    # Data-associated parameters
    diseases = ['Pso', 'AD', 'LP']  # Diseases to include
    biopsy_type = ['LESIONAL', 'NON LESIONAL']  # Biopsy_types to include
    tissue_layers = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']  # Skin layers to include in the analysis

    # Correlation analysis parameters
    radii = list(np.arange(0, 10))  # Set clustering radius
    corr_method = 'spearman'  # Correlation coefficient: spearman or pearson

    # Geneset parameters
    # Which response_types should to be tested (Th1: 'IFNG', Th2: 'IL13', Th3:'IL17');
    # can be None if response type is not needed, please check functions if they require one
    response_type = 'IL17'
    # Get genes of interest and their potential close vicinity genes
    conditional_genes, conditionalgenes_responders = gene_lists.get_publication_spatial_transcriptomics_cyto_resps(
        respone_type=response_type)

    # Path-related parameters
    task = 'Cytokines_vs_Responders'  # Name of task or project
    # Data directory inside docker: os.environ.get("DATAPATH", "/data")
    # Data directory outside docker: 'path/to/data' # e.g.'/Volumes/T7/data/annData_objects/spatial/2022-04-08'
    path_to_data = os.environ.get("DATAPATH", "/data")
    # Output directory when using docker
    path_to_output = os.getcwd()

    # Output directory
    save_folder = os.path.join(path_to_output, "results", task, str(date.today()))
    os.makedirs(save_folder, exist_ok=True)

    # Load data
    unpp_st_adata = sc.read(
        os.path.join(path_to_data, "Spatial Transcriptomics_unpp_cleaned_LPADPso_reduced.h5"))
    # Subset to specific disease(s)
    unpp_st_adata = unpp_st_adata[unpp_st_adata.obs['DISEASE'].isin(diseases)]
    # use lesion and  samples
    unpp_st_adata = unpp_st_adata[unpp_st_adata.obs['biopsy_type'].isin(biopsy_type)]

    counts_dict = dict()
    dict_corr_genes = {}

    for radius in radii:
        save_folder_tmp = os.path.join(save_folder, str(radius))
        os.makedirs(save_folder_tmp, exist_ok=True)

        print('========= Initiate Object radius = {} ========='.format(radius))
        clusterer = modul_density_clustering.SpatialDensityCluster(
            adata=unpp_st_adata, tissue_types=tissue_layers, radius=radius, corr_method=corr_method,
            conditional_genes=conditional_genes,
            conditionalgenes_responders=conditionalgenes_responders)

        print('========= Run ST Density Clustering radius = {} =========\n'.format(radius))
        df_counts, df_excluded_spot_counts, df_counts_responders, df_included_spot_counts = (
            clusterer.get_cluster_counts())

        counts_dict[radius] = df_counts
        dict_corr_genes[radius] = {}

        print('========= Calculate correlation and p-value for radius = {} =========\n'.format(radius))
        for gene in conditional_genes:
            dict_corr_genes[radius][gene] = {'pearson': [], 'spearman': []}

            proximity_genes_name = "{}_{}".format(gene, 'responder')
            weight_name = "{}_{}".format('weighted', gene)

            # Read out counts
            temp_df = df_counts[[gene, proximity_genes_name, weight_name]].copy()
            # Drop NaN values
            df_wo_nan = temp_df.dropna()

            dict_corr_genes[radius][gene] = corr_statistics.get_correlation_stats(
                df=df_wo_nan, proximity_genes_name=proximity_genes_name, gene=gene, weight_name=weight_name,
                dict_corr=dict_corr_genes[radius][gene])

        print('========= Plot correlation radius = {} =========\n'.format(radius))
        # Weighted by transcripts
        plot_spatial_correlation.plot__stwc_tissuelayers(
            df_counts=df_counts, cytokine_responders=conditionalgenes_responders,
            save_folder=save_folder_tmp, distance=radius, corr_method=corr_method,
            dict_corr=dict_corr_genes[radius])
        # Paper Figure 1G-style (all black, fixed size)
        plot_spatial_correlation.plot__stwc_tissuelayers(
            df_counts=df_counts, cytokine_responders=conditionalgenes_responders, save_folder=save_folder_tmp,
            distance=radius, corr_method=corr_method, dict_corr=dict_corr_genes[radius],
            color_dots=False, scale_dot_size=False)
        # Highlight diseases
        plot_spatial_correlation.plot__stwc_disease(
            df_counts=df_counts, cytokine_responders=conditionalgenes_responders,
            save_folder=save_folder_tmp, distance=radius, corr_method=corr_method,
            dict_corr=dict_corr_genes[radius])
        # Highlight patients
        plot_spatial_correlation.plot__stwc_patients(
            df_counts=df_counts, cytokine_responders=conditionalgenes_responders,
            save_folder=save_folder_tmp, distance=radius, corr_method=corr_method,
            dict_corr=dict_corr_genes[radius])
        # Highlight biopsy type (LESIONAL and NON LESIONAL)
        plot_spatial_correlation.plot__stwc_biopsytype(
            df_counts=df_counts, cytokine_responders=conditionalgenes_responders,
            save_folder=save_folder_tmp, distance=radius, corr_method=corr_method,
            dict_corr=dict_corr_genes[radius])

    # Read out radius, correlation, and p-value
    df_radius_vs_correlation = evaluations.build_radius_vs_correlation_df(
        dict_corr_genes=dict_corr_genes, conditional_genes=conditional_genes, response_type=response_type,
        corr_method=corr_method)
    df_respcounts_radius = evaluations.compute_responder_counts_by_radius_df(
        counts_dict=counts_dict, conditionalgenes_responders=conditionalgenes_responders, radii=radii)

    if len(radii) > 1:
        # 6. Evaluate distance via elbow plot
        evaluations.plot_evaluate_distance(
            df_spearman=df_radius_vs_correlation, cytokines=conditional_genes, min_radius=min(radii),
            save_folder=save_folder, corr_method=corr_method)

        # 7. Plot Responder counts normed by number of spots in a cluster
        evaluations.plot_responder_vs_radius(
            df_respcounts_radius=df_respcounts_radius, save_folder=save_folder)

    print('========= Save results for {} and response type {} =========\n'.format(diseases, response_type))
    # Save everything to an .xlsx file
    writer = pd.ExcelWriter(os.path.join(save_folder, 'Correlation.xlsx'), engine='xlsxwriter')
    for radius in radii:
        counts_dict[radius].to_excel(writer, sheet_name='counts_radius_{}'.format(radius))
    df_radius_vs_correlation.to_excel(writer, sheet_name='radius_vs_correlation')
    df_respcounts_radius.to_excel(writer, sheet_name='normalised_responder_counts')
    writer.close()

    print("Finished")


if __name__ == '__main__':
    main()
