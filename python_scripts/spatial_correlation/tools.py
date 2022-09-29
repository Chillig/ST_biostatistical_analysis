"""Prepare, subset, and add metaData to annData object
    File name: utils.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

import numpy as np
import anndata


def convert_categories_cytokines_responders_others(adata: anndata, cyto_responder_genes: dict):
    """Add observable column to adata marking the spots either as "Other" (2), "Responders" (1) or cytokine (0)
        - others label = 2
        - responder label = 1
        - cytokine label 0

    Parameters
    ----------
    adata : annData
    cyto_responder_genes : dict

    Returns
    -------

    """

    for cyto in cyto_responder_genes.keys():
        obs_counts = "_".join([cyto, 'counts'])
        obs_label = "_".join([cyto, 'responders'])
        obs_clusters = "_".join([cyto, 'clusters'])

        # create array with others
        if "counts" in adata.layers.keys():
            default_counts = np.copy(adata.layers['counts']).sum(axis=1)
        else:
            default_counts = np.copy(adata.X).sum(axis=1)
        default_label = np.array(['Others'] * adata.shape[0], dtype='<U64')
        default_clusters = np.array([2] * adata.shape[0])

        # Get available responder genes
        available_responders = list(set(adata.var.index) & set(cyto_responder_genes[cyto]))

        # 2. Get mask for responder genes spot positions
        # 2.1 Find index of responder genes
        varindex_responder_genes = np.where(adata.var.index[np.newaxis, :] ==
                                            np.array(available_responders)[:, np.newaxis])[1]
        # 2.2 Get count matrix of responder genes
        if "counts" in adata.layers.keys():
            responder_matrix = np.copy(adata.layers['counts'])[:, varindex_responder_genes]
        else:
            responder_matrix = np.copy(adata.X)[:, varindex_responder_genes]
        # Identify where counts > 0
        m_responders = responder_matrix > 0
        m_array_responders = m_responders.sum(axis=1).astype(bool)

        # Add label of responder
        default_label[m_array_responders] = 'Responders'
        default_clusters[m_array_responders] = 1
        default_counts[m_array_responders] = np.sum(responder_matrix, axis=1)[m_array_responders]

        # 3. Get cytokine spot positions
        if cyto in adata.var_names:
            if "counts" in adata.layers.keys():
                cyto_matrix = adata.layers['counts'][:, np.where(adata.var.index == cyto)[0]]
            else:
                cyto_matrix = adata.X[:, np.where(adata.var.index == cyto)[0]]
            m_cyto = cyto_matrix[:, 0] > 0

            # Add label cytokine positive spots
            default_clusters[m_cyto] = 0
            default_label[m_cyto] = cyto
            default_counts[m_cyto] = cyto_matrix[:, 0][m_cyto]

        # 4. Add observables to adata
        adata.obs[obs_clusters] = default_clusters
        adata.obs[obs_label] = default_label
        adata.obs[obs_counts] = default_counts

        # 5. convert to categorical
        adata.obs[obs_clusters] = adata.obs[obs_label].astype('category')
        adata.obs[obs_label] = adata.obs[obs_label].astype('category')

    return adata


def add_columns_genes(adata: anndata, genes: str, label: str, count_threshold: int = 1):
    """Add columns of counts, label, and clusters to adata object beginning with term label
        gene(s) label = 0
        others label = 1

    Parameters
    ----------
    adata : annData
        object which contains all information about the experiment
    genes : str, list
        it is possible to only provide a gene name or a list of gene names
    label : str
        label for new column
    count_threshold : int
        threshold for UMI-count: a genes need to have more than count_threshold counts to be signed with the labels

    Returns
    -------

    """
    # count threshold need to be at least 1, otherwise gene was not measured
    if count_threshold < 1:
        count_threshold = 1

    # name columns
    obs_counts = "_".join([label, 'counts'])
    obs_label = "_".join([label, 'label'])
    obs_clusters = "_".join([label, 'clusters'])

    # set default counts to 0
    default_counts = np.array([0] * adata.shape[0])
    # set default label "Others"
    default_label = np.array(['Others'] * adata.shape[0], dtype='<U64')
    default_clusters = np.array([1] * adata.shape[0])

    # Work around for permutation names of genes
    if "responder" in genes:
        genes = genes.split("_")[0]

    # Check if genes is of type list or string
    if isinstance(genes, list):
        # Get available genes in adata
        available_genes = list(set(adata.var.index) & set(genes))

        # 2. Get bool mask for gene spot positions
        # 2.1 Find index of genes
        varindex_genes = np.where(adata.var.index[np.newaxis, :] == np.array(available_genes)[:, np.newaxis])[1]
        # 2.2 Get count matrix of genes
        if "counts" in adata.layers.keys():
            genes_matrix = np.copy(adata.layers['counts'])[:, varindex_genes]
        else:
            genes_matrix = np.copy(adata.X)[:, varindex_genes]
        # Identify where counts > 0
        m_genes = genes_matrix >= count_threshold
        m_array_genes = m_genes.sum(axis=1).astype(bool)

        # Add label of gene
        default_label[m_array_genes] = label
        default_clusters[m_array_genes] = 0
        default_counts[m_array_genes] = np.sum(genes_matrix, axis=1)[m_array_genes]
    else:
        # 2. Check if gene is in adata object
        if genes in adata.var_names:
            # 2.1 Get count matrix of gene
            if "counts" in adata.layers.keys():
                gene_matrix = np.copy(adata.layers['counts'])[:, np.where(adata.var.index == genes)[0]]
            else:
                gene_matrix = np.copy(adata.X)[:, np.where(adata.var.index == genes)[0]]
            # Identify where counts > 0
            m_gene = gene_matrix[:, 0] >= count_threshold

            # 2.2 Get gene spot position
            # Add label to spots
            default_label[m_gene] = label
            default_clusters[m_gene] = 0
            default_counts[m_gene] = gene_matrix[:, 0][m_gene]

    # 3. Add observables to adata
    adata.obs[obs_clusters] = default_clusters
    adata.obs[obs_label] = default_label
    adata.obs[obs_counts] = default_counts

    # 3.3 convert to categorical
    adata.obs[obs_clusters] = adata.obs[obs_label].astype('category')
    adata.obs[obs_label] = adata.obs[obs_label].astype('category')

    return adata


def mark_spotsincluster(adata: anndata, sub_adata: anndata, spot_indices: [int, list], obs_conditional_gene_counts: str,
                        gene: str, radius: [int, str]):
    """Marking spots either as "Other" (0), "cyto- nn spots" (2) or "cyto+ spots" (1)
        - others label = 0
        - responder label = 2
        - cytokine label 1

    Parameters
    ----------
    adata : annData
    sub_adata : annData
    spot_indices: int, list
    obs_conditional_gene_counts: str
    gene : str
    radius: int, str

    Returns
    -------

    """
    # Add number of spots to index of sub_adata to read out the correct spots from adata
    # get and add index_num_spots from previous samples
    temp_sample = np.unique(sub_adata.obs['specimen'])[0]

    # check if you read out correct spots with
    ind_sample = np.where(adata.obs['specimen'] == temp_sample)[0]
    spot_indices_sample = ind_sample[spot_indices]
    if not np.all(
            sub_adata.obs.iloc[spot_indices].index == adata.obs.iloc[spot_indices_sample].index):
        print("You do not read out the correct spots!! ", np.unique(sub_adata.obs['specimen']))

    adata.obs['{}_in_sdcc_r{}'.format(gene, radius)].iloc[spot_indices_sample] = 2
    # conditional gene of interest -> label = 1
    m_cg = adata.obs[obs_conditional_gene_counts].iloc[spot_indices_sample].values > 0
    if np.any(m_cg):
        adata.obs['{}_in_sdcc_r{}'.format(gene, radius)].iloc[spot_indices_sample[m_cg]] = 1

    return adata
