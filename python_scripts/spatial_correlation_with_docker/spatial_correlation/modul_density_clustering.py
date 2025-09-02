import os

import anndata
import numpy as np
import pandas as pd
import itertools
from collections import OrderedDict
from datetime import date
import scanpy as sc

# import python scripts
import utils
import condition_graphs as cyto_graphs
import subset_data, gene_lists
import evaluations, plot_spatial_correlation
import corr_statistics as corr_stats


class ClusterDataCollector:
    def __init__(self, conditional_gene_responders: dict, responder_genes: list):
        self.cluster_rows = []
        self.responders_rows = []
        self.included_counts = []
        self.excluded_counts = []
        self.conditional_gene_responders = conditional_gene_responders
        self.responder_genes = responder_genes

    def record_cluster(self, cluster_obs: pd.DataFrame, gene: str, counts_data: dict, cluster_indices: list):
        row = {
            gene: cluster_obs[f"{gene}_counts"].sum(),
            f"{gene}_responder": cluster_obs[f"{gene}_Responders_counts"].sum(),
            "{}_{}".format('weighted', gene): cluster_obs[f"{gene}_counts"].sum(),
            "Cluster_size": (cluster_obs[f"{gene}_counts"] > 0).sum(),
            "Cluster_num_spots": len(cluster_indices),
            "disease": cluster_obs["DISEASE"].iloc[0],
            "Specimen": cluster_obs["specimen"].iloc[0],
            "biopsy_type": cluster_obs["biopsy_type"].iloc[0],
            "tissue_layer": "_".join(np.unique(cluster_obs["tissue_layer"])),
            "Patient": cluster_obs["patient"].iloc[0],
        }
        for gene in self.responder_genes:
            if gene in counts_data:
                row[gene] = counts_data[gene][cluster_indices].sum()
        self.cluster_rows.append(row)

    def add_included_counts(self, values, colname):
        self.included_counts.append((values, colname))

    def add_excluded_counts(self, values, colname):
        self.excluded_counts.append((values, colname))

    def to_dataframes(self):
        df_clusters = pd.DataFrame.from_records(self.cluster_rows)
        df_inc = pd.concat([pd.DataFrame(v, columns=[col]) for v, col in self.included_counts], ignore_index=True)
        df_exc = pd.concat([pd.DataFrame(v, columns=[col]) for v, col in self.excluded_counts], ignore_index=True)
        return df_clusters, df_exc, df_inc


class SpatialDensityCluster:
    def __init__(self, adata: anndata.AnnData, tissue_types: list,
                 conditional_genes: list, conditionalgenes_responders: dict, radius: int,
                 corr_method: str = 'pearson', find_associated_genes: bool = False, get_plots: bool = False):

        self.conditional_genes = conditional_genes
        self.conditionalgenes_responders = conditionalgenes_responders
        self.radius = radius
        self.corr_method = corr_method
        self.find_associated_genes = find_associated_genes
        self.get_plots = get_plots

        self.adata = adata
        self._prepare_adata(tissue_types=tissue_types)

    def _prepare_adata(self, tissue_types: list):
        self.adata = subset_data.data_preparation(
            self.adata, conditional_genes=self.conditional_genes,
            conditionalgenes_responders=self.conditionalgenes_responders)

        if tissue_types:
            bool_col = self.adata.obs[tissue_types] == 1
            merged = bool_col.sum(axis=1)
            self.adata = self.adata[merged == 1]

        for gen in self.conditionalgenes_responders.keys():
            self.adata.obs['{}_in_sdcc_r{}'.format(gen, self.radius)] = 0

    @staticmethod
    def index_spots(spot_coordinate: np.ndarray, d: int):
        ind_list = []
        for k in range(-d, d + 1):
            m_max = 2 * d - abs(k)
            for m in range(-m_max, m_max + 1, 2):
                if m != 0 or k != 0:
                    ind_list.append((spot_coordinate[0] + m, spot_coordinate[1] + k))
        return ind_list

    def get_location_spots(self, center_coord: np.ndarray, radius: int, all_coordinates: list):
        nn_indices = self.index_spots(center_coord, d=int(radius))
        center_idx = np.where(np.all(np.array(all_coordinates) == np.array(center_coord), axis=1))[0]
        if radius == 0:
            nn_idx = center_idx
        else:
            nn_idx = np.where(
                np.all(np.array(all_coordinates)[None, :] == np.array(nn_indices)[:, None], axis=2))[1]
        return center_idx, nn_idx

    def subgraph_nncytospots(self, subgraph, cyto_coords: list, all_coords: list, radius: int = 1):
        gene_indices, cyto_indices = [], []
        for spot in subgraph:
            cyto_idx, resp_idx = self.get_location_spots(
                center_coord=cyto_coords[spot], radius=radius, all_coordinates=all_coords)
            gene_indices.extend(np.concatenate([resp_idx, cyto_idx]))
            cyto_indices.extend(cyto_idx)
        return list(OrderedDict.fromkeys(cyto_indices)), list(OrderedDict.fromkeys(gene_indices))

    @staticmethod
    def get_responder_cytoposneg_counts(
            adata_obs: pd.DataFrame, observable_cyto, observable_resps, df, geneoffocus: str):
        col_respcytpos = f'{geneoffocus}+ responder+ spot'
        col_respcytneg = f'{geneoffocus}- responder+ spot'

        m_resp_cytopos = (adata_obs[observable_cyto] > 0) & (adata_obs[observable_resps] > 0)
        m_resp_cytoneg = (adata_obs[observable_cyto] == 0) & (adata_obs[observable_resps] > 0)

        df_temp_pos = pd.DataFrame(adata_obs[observable_resps][m_resp_cytopos].values, columns=[col_respcytpos])
        df_temp_neg = pd.DataFrame(adata_obs[observable_resps][m_resp_cytoneg].values, columns=[col_respcytneg])

        df = pd.concat([df, df_temp_pos, df_temp_neg], ignore_index=True)
        return df

    def get_cluster_counts(self):

        responder_genes = list(itertools.chain(*self.conditionalgenes_responders.values()))
        collector = ClusterDataCollector(self.conditionalgenes_responders, responder_genes)

        df_counts_responders = pd.DataFrame()  # optionally fill if needed

        specimen_names = np.unique(self.adata.obs['specimen'])
        for specimen in specimen_names:
            sub_adata = self.adata[self.adata.obs['specimen'] == specimen]
            for conditional_gene in self.conditionalgenes_responders:
                obs_cyto_counts = f"{conditional_gene}_counts"
                obs_resp_counts = f"{conditional_gene}_Responders_counts"
                obs_cyto_label = f"{conditional_gene}_label"

                cyto_spots = sub_adata[sub_adata.obs[obs_cyto_label] == conditional_gene]
                cyto_coordinates = list(zip(cyto_spots.obs['array_col'], cyto_spots.obs['array_row']))
                all_coordinates = list(zip(sub_adata.obs['array_col'], sub_adata.obs['array_row']))

                if not cyto_coordinates:
                    mask_resp = sub_adata.obs[obs_resp_counts] > 0
                    collector.add_excluded_counts(
                        sub_adata.obs[obs_resp_counts][mask_resp].values, f"{conditional_gene}_responder")
                    continue

                spot_connect = cyto_graphs.SpotConnectivity(np.array(cyto_coordinates))
                # Using KDTree: connected components
                subgraph_cyto_spots = spot_connect.get_connections(distance=2.0, tree_type="KDTree")

                graphs_cluster = []
                for subgraph in subgraph_cyto_spots:
                    cyto_idx, gene_idx = self.subgraph_nncytospots(
                        subgraph=subgraph, cyto_coords=cyto_coordinates, all_coords=all_coordinates, radius=self.radius)
                    graphs_cluster.append(gene_idx + cyto_idx)

                final_clusters = sorted(cyto_graphs.connected_components(
                    spot_connect.to_graph(graphs_cluster)), key=len, reverse=True)
                included_spots = list(itertools.chain(*final_clusters))

                counts_data = {}
                if "counts" in sub_adata.layers:
                    for gene in responder_genes:
                        if gene in sub_adata.var_names:
                            counts_data[gene] = sub_adata.layers["counts"][:, sub_adata.var_names.get_loc(gene)]
                else:
                    df = sub_adata.to_df()
                    for gene in responder_genes:
                        if gene in df.columns:
                            counts_data[gene] = df[gene].values

                for cluster in final_clusters:
                    cluster = list(cluster)
                    cluster_obs = sub_adata.obs.iloc[cluster]
                    collector.record_cluster(cluster_obs=cluster_obs, gene=conditional_gene, counts_data=counts_data,
                                             cluster_indices=cluster)

                    df_counts_responders = self.get_responder_cytoposneg_counts(
                        adata_obs=cluster_obs, observable_cyto=obs_cyto_counts, observable_resps=obs_resp_counts,
                        df=df_counts_responders, geneoffocus=conditional_gene)

                self.adata = utils.mark_spotsincluster(
                    self.adata, sub_adata, included_spots, obs_cyto_counts, conditional_gene, self.radius)

                full_index = np.arange(sub_adata.shape[0])
                mask_included = np.isin(full_index, included_spots)
                included_obs = sub_adata.obs.iloc[mask_included]
                excluded_obs = sub_adata.obs.iloc[~mask_included]

                mask_inc = included_obs[obs_resp_counts] > 0
                collector.add_included_counts(
                    included_obs[obs_resp_counts][mask_inc].values, f"{conditional_gene}_responder")

                mask_exc = excluded_obs[obs_resp_counts] > 0
                collector.add_excluded_counts(
                    excluded_obs[obs_resp_counts][mask_exc].values, f"{conditional_gene}_responder")

        df_counts, df_excluded_spot_counts, df_included_spot_counts = collector.to_dataframes()

        return df_counts, df_excluded_spot_counts, df_counts_responders, df_included_spot_counts


def example():
    today = date.today()

    # Load data:
    # Inside docker: os.environ.get("DATAPATH", "/data")
    # Outside docker: 'path/to/data' # e.g.'/Volumes/T7/data/annData_objects/spatial/2022-04-08'
    path_to_data = os.environ.get("DATAPATH", "/data")
    unpp_st_adata = sc.read(
        os.path.join(path_to_data, "Spatial Transcriptomics_unpp_cleaned_LPADPso_reduced.h5"))

    # Parameters
    diseases = ['Pso']
    response_types = ['IL17']
    tissue_layers = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']
    biopsy_type = ['LESIONAL']
    radii = [0, 1]
    corr_method = 'spearman'
    project_name = 'Cytokines_vs_Responders'

    for disease in diseases:
        disease = [disease]
        for response_type in response_types:
            # Subset to specific disease(s) and lesion samples
            mask = unpp_st_adata.obs['DISEASE'].isin(disease) & (unpp_st_adata.obs['biopsy_type'].isin(biopsy_type))
            # Subset to specific disease(s)
            unpp_st_adata_tmp = unpp_st_adata[mask]

            # store spatial information
            unpp_st_adata_tmp.obs['sample'] = unpp_st_adata_tmp.obs['sample'].cat.remove_unused_categories()

            # 1. Get cytokines and responders
            conditional_genes, conditionalgenes_responders = gene_lists.get_publication_spatial_transcriptomics_cyto_resps(
                respone_type=response_type)

            # Create Task name
            task = '{}__{}_{}'.format(project_name, "_".join(disease), response_type)

            # Create output directory
            save_folder = os.path.join(os.getcwd(), "results", task, str(today))
            os.makedirs(save_folder, exist_ok=True)

            counts_dict = dict()
            dict_corr_genes = {}
            for radius in radii:
                save_folder_tmp = os.path.join(save_folder, str(radius))
                os.makedirs(save_folder_tmp, exist_ok=True)

                print('========= Initiate Object radius = {} ========='.format(radius))
                clusterer = SpatialDensityCluster(
                    adata=unpp_st_adata_tmp, tissue_types=tissue_layers, radius=radius, corr_method=corr_method,
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

                    dict_corr_genes[radius][gene] = corr_stats.get_correlation_stats(
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

            print('========= Save results for {} and response type {} =========\n'.format(disease, response_type))
            # Save everything to an .xlsx file
            writer = pd.ExcelWriter(os.path.join(save_folder, 'Correlation.xlsx'), engine='xlsxwriter')
            for radius in radii:
                counts_dict[radius].to_excel(writer, sheet_name='counts_radius_{}'.format(radius))
            df_radius_vs_correlation.to_excel(writer, sheet_name='radius_vs_correlation')
            df_respcounts_radius.to_excel(writer, sheet_name='normalised_responder_counts')
            writer.close()

    print("Finished")


if __name__ == '__main__':
    example()
