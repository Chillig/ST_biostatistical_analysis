"""Calculate (Weighted) Spatial Pearson Correlation for different methods based on conditional nearest neighbors
    File name: density_clustering.py
    Author: Christina Hillig, Ali Farnoud-Niedermayr
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

# import scripts
from python_scripts.spatial_correlation import tools, condition_graphs as cyto_graphs, helper_functions as ht, \
    refine_responder_genes
from python_scripts.spatial_correlation.plots import plot_clusters, plot_count_distributions, plot_evaluations, \
    plot_spatial_correlation
from python_scripts.spatial_correlation import helper_functions

import scanpy as sc
import numpy as np
import pandas as pd
import os
from collections import OrderedDict
import itertools


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
    """Get index of spots in all_coordinates from center spot and nn spots

    Parameters
    ----------
    center_coord : tuple
        coordinate of center spot
    distance : int
        radius starting from center spot
    all_coordinates : list of tuples
        list containing all coordinates

    Returns
    -------

    """
    # 1. Get index of nn spots
    nn_indices = index_spots(center_coord, d=int(distance))

    # 2. Get index of nn spots in list of all coordinates and of cyto+ spots
    # 2.1 Index of cyto+ spots
    varindex_center_gene = np.where(np.all(np.array(all_coordinates) == np.array(center_coord), axis=1))[0]
    # 2.2 Indices of nn responder spots
    if distance == 0:
        varindex_nn_genes = varindex_center_gene
    else:
        varindex_nn_genes = np.where(
            np.all(np.array(all_coordinates)[np.newaxis, :] == np.array(nn_indices)[:, np.newaxis], axis=2))[1]

    # Check if above command extracts nn spots positions in adata object -> yep, it does :)
    # np.array(all_coordinates)[varindex_responder_genes] == np.array(nn_indices)nn_indices

    return varindex_center_gene, varindex_nn_genes


def subgraph_nncytospots(subgraph, cyto_coordinates, all_coordinates, distance=1):
    """Read out cytokine and responder spots per sub-graph

    Parameters
    ----------
    subgraph : set
        nodes (cyto+ spots) in subgraph
    cyto_coordinates : list of tuples
        cytokine coordinates stored as tuples in list [(x1, y1), ...]
    all_coordinates : list of tuples
        all coordinates stored as tuples in list [(x1, y1), ...]
    distance : int
        radius

    Returns
    -------

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
    """ Read out responder counts in cyto+ spots and cyto- spots in the clusters
    -> find evidence that counts of responder are highest in cyto+ spots (boxplot + stats test)

    Parameters
    ----------
    adata_obs : pd.DataFrame
        contains all columns of adata of specific cluster
    obervable_cyto : str
        name of observable to read out cytokine counts
    obervable_resps : str
        name of observable to read out responder counts
    df : pd.DataFrame
        store information of responder counts in cyto+ and cyto- spots in a specific cluster
    cytokine : str
        name of cytokine

    Returns
    -------

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


def get_cluster_counts(adata, tissue_types, cytokine_responders, save_folder, distance=1, find_responders=False,
                       get_plots=False):
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
    find_responders : bool
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
        # Add obs to adata stating whether a spot belongs to a cluster or not
        # init with assumption: no spots are in a cluster
        adata.obs['{}_in_sdcc_r{}'.format(gen, distance)] = 0

        col_names_resps.extend(["".join([gen, '+ responder+ spot']), "".join([gen, '- responder+ spot'])])
        # save colnames for df_counts
        col_names.append("_".join([gen, 'responder']))
    col_names.extend(['disease', 'Cluster_size', 'Cluster_num_spots', 'Specimen'])
    # For all genes - should be done only for normed counts ..
    # col_names.extend(bg_genes)
    # Add name of responder genes
    responder_genes = list(itertools.chain(*list(cytokine_responders.values())))
    col_names.extend(responder_genes)
    col_names = list(np.unique(col_names))
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
    specimen_names = np.unique(adata.obs['specimen'])
    # Get correlation per sample -- for loop
    index_counter = 0
    for ind_specimen, specimen in enumerate(specimen_names):
        # get adata of sample
        sub_adata = adata[adata.obs['specimen'] == specimen].copy()
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
                index_counter_specimen = 0
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
                    # Read out counts for each responder gene individually
                    for resp_gene in responder_genes:
                        if resp_gene in sub_adata.var_names:
                            if find_responders and ("counts" in sub_adata.layers.keys()):
                                counts_respgene = np.copy(
                                    sub_adata.layers["counts"])[:, np.where(sub_adata.var.index == resp_gene)[0]]
                                # create mask
                                m_respgene_cytopos = counts_respgene > 0
                                df_counts.at[index_counter, resp_gene] = np.divide(
                                    counts_respgene[m_respgene_cytopos].sum(), 1)

                            elif "counts" in sub_adata.layers.keys():
                                counts_respgene = np.copy(adata.X)[:, np.where(adata.var.index == resp_gene)[0]]
                                # create mask
                                m_respgene_cytopos = counts_respgene > 0
                                df_counts.at[index_counter, resp_gene] = np.divide(
                                    counts_respgene[m_respgene_cytopos].sum(), 1)

                            else:
                                counts_respgene = sub_adata.to_df()
                                # build mask to read out counts from gene+ spots
                                m_respgene_cytopos = counts_respgene[resp_gene] > 0
                                df_counts.at[index_counter, resp_gene] = np.divide(
                                    counts_respgene.iloc[list(indices_spots), :][resp_gene][
                                        m_respgene_cytopos].sum(), 1)
                        else:
                            print("Responder gene not in dataset: ", resp_gene)

                    # np.count_nonzero(nn_obs[obs_resp_counts])
                    df_counts.at[index_counter, 'disease'] = np.unique(sub_adata.obs['DISEASE'])[0]
                    df_counts.at[index_counter, 'Cluster_size'] = np.count_nonzero(cluster_obs[obs_cyto_counts])
                    # This assumes that each cyto+ spot contains any of the responder genes
                    # should be equal to len(list(indices_spots)) == np.count_nonzero(cluster_obs[obs_resp_counts])
                    # BUT its not as some responder genes are not found in cyto+ spot ..
                    df_counts.at[index_counter, 'Cluster_num_spots'] = len(list(indices_spots))
                    df_counts.at[index_counter, 'tissue_layer'] = "_".join(np.unique(cluster_obs['tissue_layer']))
                    df_counts.at[index_counter, 'Specimen'] = specimen
                    df_counts.at[index_counter, 'Patient'] = np.unique(sub_adata.obs['patient'])[0]

                    # Read out responder counts in cyto+ spots and cyto- spots in the clusters
                    df_counts_responders = get_responder_cytoposneg_counts(
                        adata_obs=cluster_obs, obervable_cyto=obs_cyto_counts, obervable_resps=obs_resp_counts,
                        df=df_counts_responders, cytokine=cyto)

                    index_counter = index_counter + 1
                    index_counter_specimen = index_counter_specimen + 1

                # Save all genes in close proximity of cyto+ spots
                adata = tools.mark_spotsincluster(
                    adata=adata, sub_adata=sub_adata, spot_indices=included_spots,
                    obs_conditional_gene_counts=obs_cyto_counts, gene=cyto, radius=distance)

                if get_plots:
                    plot_count_distributions.plot_counts(
                        dfx=df_counts["_".join([cyto, 'responder'])], dfy=df_counts[cyto],
                        index_counter=index_counter, index_counter_sample=index_counter_specimen, cyto=cyto,
                        sample=specimen, save_folder=save_folder, distance=distance)

                    # Plot on image clusters of cyto+ and nn responder spots
                    sub_dict = {k: cytokine_responders[k] for k in cytokine_responders.keys() & {cyto}}
                    sevenspots_adata = sub_adata[included_spots]
                    sevenspots_adata = tools.convert_categories_cytokines_responders_others(
                        adata=sevenspots_adata, cyto_responder_genes=sub_dict)
                    obs_label = "_".join([cyto, 'responders'])
                    obs_counts = "_".join([cyto, 'counts'])
                    samples, crops_img = helper_functions.get_cropped_sampleimg(img_key='highres')
                    if sevenspots_adata.obs['sample'].cat.categories in samples:
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

    return df_counts, df_excluded_spot_counts, df_counts_responders, df_included_spot_counts, adata


def run_spatialcorrelation(adata, tissue_types, cytokine_responders, save_folder, radius, sig, corr_method,
                           find_responders=False, get_plots=False):
    """Perform spatial Correlation Analysis and create Plots

    Parameters
    ----------
    adata : annData
    tissue_types : list
    cytokine_responders : dict
    save_folder : str
    radius : int
    sig : None, list
    corr_method : str
    find_responders : bool
    get_plots : bool

    Returns
    -------

    """

    # Create new folder in this directory for specific distance
    save_folder = os.path.join(save_folder, str(radius))
    os.makedirs(save_folder, exist_ok=True)

    # Get Density clusters and counts of cytokines and responders
    df_counts, df_excluded_spot_counts, df_counts_responders, df_included_spots, adata = get_cluster_counts(
        adata=adata, tissue_types=tissue_types, cytokine_responders=cytokine_responders,
        distance=radius, save_folder=save_folder, find_responders=find_responders, get_plots=get_plots)
    # Info: BET3L and FAM26D are no longer in dataset.. -> overall 71 responder genes

    # Perform DGE analysis between cyto and nn spots
    if find_responders and ("counts" in adata.layers.keys()):
        refine_responder_genes.rank_cyto_vs_others_genes(adata, cytokine_responders, save_folder, radius=radius)

    # 4. Plot correlation
    # Weighted by transcripts
    figs4_wtcorr_pval = plot_spatial_correlation.plot__stwc_tissuelayers(
        df_counts=df_counts, cytokine_responders=cytokine_responders, save_folder=save_folder, distance=radius,
        corr_method=corr_method)
    sig.append(figs4_wtcorr_pval)
    _ = plot_spatial_correlation.plot__stwc_patients(
        df_counts=df_counts, cytokine_responders=cytokine_responders, save_folder=save_folder,
        distance=radius, corr_method=corr_method)

    # Plot Workflow figure 1G
    if radius == 1:
        plot_spatial_correlation.plot__stwc_tissuelayers_figure1g(
            df_counts=df_counts, save_folder=save_folder, distance=radius, corr_method=corr_method)

    return sig, df_counts, adata


def data_preparation(adata, tissue_types, epidermis_layers, conditional_genes, conditionalgenes_responders):
    """Prepare data before applying conditional density clustering algorithm

    Parameters
    ----------
    adata : annData
    tissue_types : list
    conditional_genes : list
    epidermis_layers : str, list
    conditionalgenes_responders : dict

    Returns
    -------

    """
    # Subset adata to tissue_types of interest: upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS'
    if tissue_types:
        bool_col = adata.obs[tissue_types] == 1
        merged = bool_col.sum(axis=1)
        adata = adata[merged == 1].copy()
        # Rename tissue region 'JUNCTION' to basal EPIDERMIS because some spots got both labels
        adata = ht.interface_to_epidermis(adata, tissue_layers='tissue_layer', epidermis_layers=epidermis_layers)

    # Get counts of cyotkines and their responders in the EPIDERMIS
    # - distance between spots: 100 µm, spot diameter: 55 µm
    # - better: use index array
    resp_label = []
    max_nl_counts_cytokines = dict()
    for cyto in conditional_genes:
        adata = tools.add_columns_genes(adata=adata, genes=cyto, label=cyto, count_threshold=1)
        # Calculate new cut-off using non lesion as reference
        max_nl_counts_cytokines[cyto] = adata.obs[
            '{}_counts'.format(cyto)][adata.obs['biopsy_type'] == 'NON LESIONAL'].max()

    for ind, cyto_resps in enumerate(conditionalgenes_responders.keys()):
        resp_label.append("_".join([cyto_resps, "Responders"]))
        adata = tools.add_columns_genes(adata=adata, genes=conditionalgenes_responders[cyto_resps],
                                        label=resp_label[ind], count_threshold=1)
        # Calculate new cut-off using non lesion as reference
        max_nl_counts_cytokines[resp_label[ind]] = adata.obs[
            '{}_counts'.format(resp_label[ind])][adata.obs['biopsy_type'] == 'NON LESIONAL'].max()

    print(max_nl_counts_cytokines)

    return adata


def main(adata, save_folder, tissue_types, epidermis_layers, radii, corr_method, conditional_genes,
         conditionalgenes_responders, find_responders=False, get_plots=False):
    """

    Parameters
    ----------
    adata : annData
    save_folder : str
    tissue_types : str, list
    epidermis_layers : str, list
    radii : int, list
    corr_method : str
    conditional_genes : list
    conditionalgenes_responders : dict
    find_responders: bool
    get_plots : bool

    Returns
    -------

    """

    # 1. prepare adata object
    adata = data_preparation(adata=adata, tissue_types=tissue_types, conditional_genes=conditional_genes,
                             conditionalgenes_responders=conditionalgenes_responders, epidermis_layers=epidermis_layers)

    # 3. Run conditional clustering and calculate (spatially weighted) correlation
    sig = []
    counts_dict = dict()
    if isinstance(radii, list):
        for radius in radii:
            sig, df_counts, adata = run_spatialcorrelation(
                adata=adata, tissue_types=tissue_types, cytokine_responders=conditionalgenes_responders,
                save_folder=save_folder, radius=radius, sig=sig, get_plots=get_plots, corr_method=corr_method,
                find_responders=find_responders)
            counts_dict[radius] = df_counts
    else:
        sig, df_counts, adata = run_spatialcorrelation(
            adata=adata, tissue_types=tissue_types, cytokine_responders=conditionalgenes_responders,
            save_folder=save_folder, radius=radii, sig=sig, get_plots=get_plots, corr_method=corr_method,
            find_responders=find_responders)
        counts_dict[radii] = df_counts

    # Boxplot of optimal radius of responder genes in density clusters (2 & LESIONAL) vs L (1 & LESIONAL)
    df_stats_responders_in_vs_outlesion_sdc = plot_evaluations.plot_responder_in_sdc_outside_l_tissue(
        adata=adata, conditionalgenes_responders=conditionalgenes_responders, save_folder=save_folder)
    df_stats_cytokines_responders_in_sdc = plot_evaluations.plot_in_sdc_cytokine_vs_responder(
        adata=adata, conditionalgenes_responders=conditionalgenes_responders, save_folder=save_folder)

    if ('counts' in adata.layers) & (find_responders is True):
        # Save adata obj for DGE analysis to identify cytokine related genes Fig. S8
        sc.write(os.path.join(save_folder, 'SDC_adata.h5'), adata)

    if len(sig) > 1:
        # 6. Evaluate distance via elbow plot
        df_radius_vs_spearman = plot_evaluations.plot_evaluate_distance(
            significance=sig, cytokines=conditional_genes, min_radius=min(radii), save_folder=save_folder)

        # 7. Plot Responder counts normed by number of spots in a cluster
        plot_evaluations.plot_responder_vs_radius(
            counts_dict=counts_dict, conditionalgenes_responders=conditionalgenes_responders,
            radii=radii, save_folder=save_folder)
    else:
        df_radius_vs_spearman = []

    return adata, counts_dict, df_stats_responders_in_vs_outlesion_sdc, \
           df_stats_cytokines_responders_in_sdc, df_radius_vs_spearman
