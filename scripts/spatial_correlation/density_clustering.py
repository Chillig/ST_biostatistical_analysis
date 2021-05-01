"""Calculate (Weighted) Spatial Pearson Correlation for different methods based on conditional nearest neighbors
    File name: density_clustering.py
    Author: Christina Hillig, Ali Farnoud-Niedermayr
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

# import scripts
from scripts.spatial_correlation import tools, condition_graphs as cyto_graphs
from scripts.spatial_correlation import plot_clusters, plot_spatial_correlation, plot_evaluations, \
    plot_count_distributions
from scripts.utils import gene_lists

import numpy as np
import pandas as pd
import os
from collections import OrderedDict


# -------------------------------------------------------------------------------------------------------------------- #
def index_spots(spot_coordinate, d):
    """Get index of nearest neighbor spots
        k: index of row
        m: index of column

    Parameters
    ----------
    spot_coordinate : numpy.array
        coordinate of spot of interest
    d : int
        distance of spot of interest to its next neighbor spots

    Returns
    -------

    """

    ind_list = []
    for k in range(-d, d + 1):
        m_max = 2 * d - np.abs(k)
        for m in range(-m_max, m_max + 1, 2):
            if m == 0 and k == 0:
                continue
            else:
                ind_list.append((spot_coordinate[0] + m, spot_coordinate[1] + k))

    return ind_list


def get_location_spots(center_coord, distance, all_coordinates):
    """Get index of spots in all_coordiantes from center spot and nn spots

    Parameters
    ----------
    center_coord : tuple
        coordinate of center spot
    distance : int
        radius starting from center spot
    all_coordinates : list of tuples
        list conatining all coordinates

    Returns
    -------

    """
    # 1. Get index of nn spots
    nn_indices = index_spots(center_coord, d=int(distance))

    # 2. Get index of nn spots in list of all coordinates and of cyto+ spots
    # 2.1 Index of cyto+ spots
    varindex_center_gene = np.where(np.all(np.array(all_coordinates) == np.array(center_coord), axis=1))[0]
    # 2.2 Indices of nn responder spots
    varindex_nn_genes = np.where(
        np.all(np.array(all_coordinates)[np.newaxis, :] == np.array(nn_indices)[:, np.newaxis], axis=2))[1]

    # Check if above command extracts nn spots positions in adata object -> yep, it does :)
    # np.array(all_coordinates)[varindex_responder_genes] == np.array(nn_indices)nn_indices

    return varindex_center_gene, varindex_nn_genes


def subgraph_nncytospots(subgraph, cyto_coordinates, all_coordinates, distance=1):
    """
    Read out cytokine and responder spots per sub-graph

    :param subgraph: [set] nodes (cyto+ spots) in subgraph
    :param cyto_coordinates: [list of tuples] cytokine coordinates stored as tuples in list [(x1, y1), ...]
    :param all_coordinates: [list of tuples] all coordinates stored as tuples in list [(x1, y1), ...]
    :param distance: [int]
    :return:
    """
    # Get nn spots and remove duplicates
    varindex_genes = []
    varindex_cytos = []
    for cc_spot in subgraph:
        # 4.
        varindex_cyto_gene, varindex_responder_genes = get_location_spots(
            center_coord=cyto_coordinates[cc_spot], all_coordinates=all_coordinates, distance=distance)
        # 5. Read out observations of nn spots and spot of interest (cytokine+)
        # add also center spot because it probably contains also responder genes
        temp_varindex_responder_genes = np.concatenate([varindex_responder_genes,
                                                        varindex_cyto_gene])
        varindex_genes.extend(temp_varindex_responder_genes)
        varindex_cytos.extend(varindex_cyto_gene)

    # remove duplicates in list
    varindex_genes = list(OrderedDict.fromkeys(varindex_genes))

    return varindex_cytos, varindex_genes


def get_responder_cytoposneg_counts(adata_obs, obervable_cyto, obervable_resps, df, cytokine):
    """
    Read out responder counts in cyto+ spots and cyto- spots in the clusters
    -> find evidence that counts of responder are highest in cyto+ spots (boxplot + stats test)

    :param adata_obs: [dataframe] contains all columns of adata of specific cluster
    :param obervable_cyto: [string] name of observable to read out cytokine counts
    :param obervable_resps: [string] name of observable to read out responder counts
    :param df: [pandas.Dataframe] store information of responder counts in cyto+ and cyto- spots in a specific cluster
    :param cytokine: [string] name of cytokine
    :return:
    """

    # create column name
    col_respcytpos = ''.join([cytokine, "+ responder+ spot"])
    col_respcytneg = ''.join([cytokine, "- responder+ spot"])
    # mask to read out responder counts in cytokine+ spots
    m_resp_cytopos = (adata_obs[obervable_cyto] > 0) & (adata_obs[obervable_resps] > 0)
    # mask to read out only responder counts from cyto- but responder+ spots
    m_resp_cytoneg = (adata_obs[obervable_cyto] == 0) & (adata_obs[obervable_resps] > 0)

    # Add counts to data frame
    df_temp = pd.DataFrame(adata_obs[obervable_resps][m_resp_cytopos].values,
                           columns=[col_respcytpos])
    df = df.append(df_temp, ignore_index=True)
    df_temp = pd.DataFrame(adata_obs[obervable_resps][m_resp_cytoneg].values,
                           columns=[col_respcytneg])
    df = df.append(df_temp, ignore_index=True)

    return df


def get_cluster_counts(adata, tissue_types, cytokine_responders, save_folder, distance=1, get_plots=False):
    """Get counts of nearest neighbor and center spot

    Parameters
    ----------
    adata : annData
    tissue_types : list of str
        contains the manually annotated tissue layers used in the analysis
    cytokine_responders : dict
        Dictionary containing cytokines as keys and responder genes as values
    save_folder : str
        path where output shall be saved
    distance : int
        max. distance of neighbouring spots to center spot
    get_plots : bool
        whether to plot in-between-steps

    Returns
    -------

    """

    # Subset adata such that only tissue types of interest are in it
    if tissue_types:
        bool_col = adata.obs[tissue_types] == 1
        merged = bool_col.sum(axis=1)
        adata = adata[merged == 1]

    # Create Dataframe to save counts of cytokine and its responder gene
    col_names = list(cytokine_responders.keys())
    col_names_resps = []
    for gen in cytokine_responders.keys():
        col_names_resps.extend(["".join([gen, '+ responder+ spot']), "".join([gen, '- responder+ spot'])])
        # save colnames for df_counts
        col_names.append("_".join([gen, 'responder']))
    col_names.extend(['disease', 'Cluster_size'])
    df_counts = pd.DataFrame(columns=col_names)

    # create dataframe to save repsonder counts in cyto+ and cyto- spots in clusters -> can this be more efficient?
    df_counts_responders = pd.DataFrame(columns=col_names_resps)

    # Dataframe to get the distribution of responder spots which are excluded in the analysis
    col_names = []
    for gen in cytokine_responders.keys():
        col_names.append("_".join([gen, 'responder']))
    df_excluded_spot_counts = pd.DataFrame(columns=col_names)
    df_included_spot_counts = pd.DataFrame(columns=col_names)

    # Get name of all samples
    samples = np.unique(adata.obs['sample'])
    # Get correlation per sample -- for loop
    index_counter = 0
    for ind_sample, sample in enumerate(samples):
        # get adata of sample
        sub_adata = adata[adata.obs['sample'] == sample].copy()
        # per cytokine -- for loop
        for cyto in cytokine_responders.keys():
            # get colnames which contain the counts and labels oth a cytokine, its responder genes and others
            obs_cyto_counts = "_".join([cyto, 'counts'])
            obs_resp_counts = "_".join([cyto, 'Responders', 'counts'])
            obs_cyto_label = "_".join([cyto, 'label'])

            # Method 1: get nearest neighbor (nn) indices of only cytokine+ spot and find nn responder+ spots;
            # merge clusters if they share responder+ spots --> final method used in publication
            # 1. Read out all cyto+ spots
            cyto_spots = sub_adata[sub_adata.obs[obs_cyto_label] == cyto]
            # 2. Save coordinates of cyto+ spots of a sample in list of tuples [(x1, y1), (x2, y2), ...]
            cyto_coordinates = list(zip(*(cyto_spots.obs['array_col'].values, cyto_spots.obs['array_row'].values)))
            # 3. Save coordinates of ALL spots of a sample in list of tuples [(x1, y1), (x2, y2), ...]
            all_coordinates = list(zip(*(sub_adata.obs['array_col'].values, sub_adata.obs['array_row'].values)))

            if len(cyto_coordinates) > 0:
                graphs_cluster = []
                # 4.1 Get all cyto+ spots which cluster together in a graph by distance d=2
                subgraph_cyto_spots = cyto_graphs.get_connections(coordinates=np.array(cyto_coordinates),
                                                                  distance=2)

                # 4.2 loop -- read out responder and cytokine counts per sub-graph
                index_counter_sample = 0
                for ind_subgraph, subgraph in enumerate(subgraph_cyto_spots):
                    # Read out nn responder spots and add to graph
                    varindex_cytos, varindex_genes = subgraph_nncytospots(
                        subgraph=subgraph, cyto_coordinates=cyto_coordinates, all_coordinates=all_coordinates,
                        distance=distance)

                    # save indices of spots belonging to same cluster in a list, dict or whatever?
                    indices_subgraph_cluster = varindex_genes.copy()
                    indices_subgraph_cluster.extend(varindex_cytos)
                    graphs_cluster.append(indices_subgraph_cluster)

                # 4.3 Check if clusters overlap
                graph_tree_clusters = cyto_graphs.to_graph(graphs_cluster)
                final_cluster_indices = sorted(cyto_graphs.connected_components(graph_tree_clusters),
                                               key=len, reverse=True)

                # Create list of spots which are included in the correlation analysis
                included_spots = []
                # 5. Read out counts with indices
                for indices_spots in final_cluster_indices:
                    included_spots.extend(list(indices_spots))
                    cluster_obs = sub_adata.obs.iloc[list(indices_spots)]
                    # Read out counts and save in dataframe
                    # Don't normalise because otherwise the size/density of cytokine+ spots would be neglected
                    # build mask to read out counts only at those spots which are cytokine+ -- redundant step
                    m_cytos = cluster_obs[obs_cyto_counts] > 0
                    df_counts.at[index_counter, cyto] = np.divide(cluster_obs[obs_cyto_counts][m_cytos].sum(), 1)
                    # build mask to read out counts from responder+ spots -- redundant step
                    m_resp_cytopos = cluster_obs[obs_resp_counts] > 0
                    df_counts.at[index_counter, "_".join([cyto, 'responder'])] = np.divide(
                        cluster_obs[obs_resp_counts][m_resp_cytopos].sum(), 1)
                    # np.count_nonzero(nn_obs[obs_resp_counts])
                    df_counts.at[index_counter, 'disease'] = np.unique(sub_adata.obs['disease'])[0]
                    df_counts.at[index_counter, 'Cluster_size'] = np.count_nonzero(cluster_obs[obs_cyto_counts])
                    df_counts.at[index_counter, 'tissue_layer'] = "_".join(np.unique(cluster_obs['tissue_type']))

                    # Read out responder counts in cyto+ spots and cyto- spots in the clusters
                    df_counts_responders = get_responder_cytoposneg_counts(
                        adata_obs=cluster_obs, obervable_cyto=obs_cyto_counts, obervable_resps=obs_resp_counts,
                        df=df_counts_responders, cytokine=cyto)

                    index_counter = index_counter + 1
                    index_counter_sample = index_counter_sample + 1

                if get_plots:
                    plot_count_distributions.plot_counts(
                        dfx=df_counts["_".join([cyto, 'responder'])], dfy=df_counts[cyto],
                        index_counter=index_counter, index_counter_sample=index_counter_sample, cyto=cyto,
                        sample=sample, save_folder=save_folder, distance=distance)

                    # Plot on image clusters of cyto+ and nn responder spots
                    sub_dict = {k: cytokine_responders[k] for k in cytokine_responders.keys() & {cyto}}
                    sevenspots_adata = sub_adata.copy()[included_spots]
                    sevenspots_adata = tools.convert_categories_cytokines_responders_others(
                        adata=sevenspots_adata, cyto_responder_genes=sub_dict)
                    obs_label = "_".join([cyto, 'responders'])
                    obs_counts = "_".join([cyto, 'counts'])
                    plot_clusters.plot_clusters_counts(
                        sevenspots_adata, cyto, save_folder=save_folder, obs_label=obs_label, obs_counts=obs_counts)

                # 11. Investigate distribution of all responder spots which are excluded in the correlation analysis
                mask_included_spot_indices = np.in1d(np.arange(0, sub_adata.shape[0], 1), included_spots)
                excluded_adata_obs = sub_adata.obs.iloc[~mask_included_spot_indices]
                # mask to read out only responder+ spot counts
                m_resp = excluded_adata_obs[obs_resp_counts] > 0
                df_temp = pd.DataFrame(excluded_adata_obs[obs_resp_counts][m_resp].values,
                                       columns=["_".join([cyto, 'responder'])])
                df_excluded_spot_counts = df_excluded_spot_counts.append(df_temp, ignore_index=True)

                included_adata_obs = sub_adata.obs.iloc[mask_included_spot_indices]
                # mask to read out only responder+ spot counts
                m_resp = included_adata_obs[obs_resp_counts] > 0
                df_temp = pd.DataFrame(included_adata_obs[obs_resp_counts][m_resp].values,
                                       columns=["_".join([cyto, 'responder'])])
                df_included_spot_counts = df_included_spot_counts.append(df_temp, ignore_index=True)
            else:
                # Add counts of responder spots
                # 11. Investigate distribution of all responder spots which are excluded in the correlation analysis
                m_resp = sub_adata.obs[obs_resp_counts] > 0
                df_temp = pd.DataFrame(sub_adata.obs[obs_resp_counts][m_resp].values,
                                       columns=["_".join([cyto, 'responder'])])
                df_excluded_spot_counts = df_excluded_spot_counts.append(df_temp, ignore_index=True)

    return df_counts, df_excluded_spot_counts, df_counts_responders, df_included_spot_counts


def run_spatialcorrelation(adata, tissue_types, cytokine_responders, save_folder, radius, sig, get_plots=False):
    """Perform spatial Correlation Analysis and create Plots

    Parameters
    ----------
    adata : annData
    tissue_types : list
    cytokine_responders : dict
    save_folder : str
    radius : int
    sig : None, list
    get_plots : bool

    Returns
    -------

    """

    # Create new folder in this directory for specific distance
    save_folder = os.path.join(save_folder, str(radius))
    os.makedirs(save_folder, exist_ok=True)

    # Get Density clusters and counts of cytokines and responders
    df_counts, df_excluded_spot_counts, df_counts_responders, df_included_spots = get_cluster_counts(
        adata=adata, tissue_types=tissue_types, cytokine_responders=cytokine_responders,
        distance=radius, save_folder=save_folder, get_plots=get_plots)

    # 4. Plot correlation
    # Pearson Correlation -> Workflow Figure 1
    _ = plot_spatial_correlation.plot__spatial_correlation(
        df_counts=df_counts, cytokine_responders=cytokine_responders, save_folder=save_folder, distance=radius)
    # Weighted Pearson Correlation by number of cytokines in cluster -> Figure 4C
    corr_pval = plot_spatial_correlation.plot__spatial_weighted__correlation_wotissuelabels(
        df_counts=df_counts, cytokine_responders=cytokine_responders, save_folder=save_folder, distance=radius)

    # 5.1 Plot distribution of responder spot counts which a excluded in the analysis
    plot_count_distributions.plot_excluded_responder_spots(
        df_spot_counts=df_excluded_spot_counts, cytokines=list(cytokine_responders.keys()), save_folder=save_folder)
    plot_count_distributions.plot_included_responder_spots(
        df_spot_counts=df_included_spots, cytokines=cytokine_responders.keys(), save_folder=save_folder)
    # 5.2 Plot distribution of responder spot counts
    plot_count_distributions.plot_distribution_respondercounts(
        adata=adata, t_cell_cytocines=list(cytokine_responders.keys()), save_folder=save_folder)

    # 5.3 Plot distribution of responder in cyto+ spots against responder counts in cyto- spots
    plot_count_distributions.plot_responder_boxplot(
        df=df_counts_responders, cytokines=cytokine_responders.keys(), save_folder=save_folder)

    plot_count_distributions.plot_compare_inex_respondercounts(
        adata=adata, df_includedcounts=df_included_spots, df_excludedcounts=df_excluded_spot_counts,
        cytokines=cytokine_responders.keys(), save_folder=save_folder)

    plot_count_distributions.plot_responder_distributions(
        adata=adata, df_included_responder=df_included_spots, df_excluded_responder=df_excluded_spot_counts,
        t_cell_cytocines=cytokine_responders.keys(), save_folder=save_folder)

    sig.append(corr_pval)

    return sig


def data_preparation(adata, tissue_types, conditional_genes, conditionalgenes_responders):
    """Prepare data before applying conditional density clustering algorithm

    Parameters
    ----------
    adata : annData
    tissue_types : list
    conditional_genes : list
    conditionalgenes_responders : dict

    Returns
    -------

    """
    # prepare adata object
    adata = tools.prepare_rawadata(adata=adata)

    # Subset adata to tissue_types of interest: upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS'
    if tissue_types:
        bool_col = adata.obs[tissue_types] == 1
        merged = bool_col.sum(axis=1)
        adata = adata[merged == 1]
        # Rename tissue region 'INTERFACE' to basal EPIDERMIS because some spots got both labels
        m_interface = adata.obs['tissue_type'] == 'INTERFACE'
        adata.obs['tissue_type'][m_interface] = 'basal EPIDERMIS'

    # Get counts of cyotkines and their responders in the EPIDERMIS
    # - distance between spots: 100 µm, spot diameter: 55 µm
    # - better: use index array
    resp_label = []
    for cyto in conditional_genes:
        adata = tools.add_columns_genes(adata=adata, genes=cyto, label=cyto)
    for ind, cyto_resps in enumerate(conditionalgenes_responders.keys()):
        resp_label.append("_".join([cyto_resps, "Responders"]))
        adata = tools.add_columns_genes(adata=adata, genes=conditionalgenes_responders[cyto_resps],
                                        label=resp_label[ind])

    return adata


def main(adata, save_folder, tissue_types, radii, get_plots=False):
    """

    Parameters
    ----------
    adata : annData
    save_folder : str
    tissue_types : str, list
    radii : int, list
    get_plots : bool

    Returns
    -------

    """
    # 1. Get cytokines and responders
    conditional_genes, _, conditionalgenes_responders = gene_lists.get_publication_cyto_resps()

    # 2. prepare adata object
    adata = data_preparation(adata=adata, tissue_types=tissue_types, conditional_genes=conditional_genes,
                             conditionalgenes_responders=conditionalgenes_responders)

    # 3. Run conditional clustering and calculate (spatially weighted) correlation
    sig = []
    if isinstance(radii, list):
        for radius in radii:
            sig = run_spatialcorrelation(adata=adata, tissue_types=tissue_types,
                                         cytokine_responders=conditionalgenes_responders, save_folder=save_folder,
                                         radius=radius, sig=sig, get_plots=get_plots)
    else:
        sig = run_spatialcorrelation(adata=adata, tissue_types=tissue_types,
                                     cytokine_responders=conditionalgenes_responders,
                                     save_folder=save_folder, radius=radii, sig=sig, get_plots=get_plots)

    if len(sig) > 1:
        # 6. Evaluate distance via elbow plot
        plot_evaluations.plot_evaluate_distance(
            significance=sig, cytokines=conditional_genes, save_folder=save_folder)
