"""Figure parameters
    File name: helper_functions.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""
import anndata
import numpy as np
from operator import itemgetter


def remove_junction_tissuelayer_spottype(adata: anndata.AnnData) -> anndata.AnnData:
    # remove JUNCTION in tissue_layer
    adata = rename_junction_in_tissuelayers(adata=adata)

    # remove JUNCTION in spot_type
    adata = rename_junction_in_spottypes(adata=adata)

    return adata


def rename_junction_in_tissuelayers(adata: anndata.AnnData) -> anndata.AnnData:
    # remove JUNCTION
    adata.obs["tissue_layer"] = adata.obs["tissue_layer"].astype(str)
    mask = adata.obs["tissue_layer"] == "JUNCTION"
    mask_middle_epidermis = adata.obs["middle EPIDERMIS"] == 1
    mask_basal_epidermis = adata.obs["basal EPIDERMIS"] == 1
    mask_upper_epidermis = adata.obs["upper EPIDERMIS"] == 1
    mask_dermis_1 = adata.obs["DERdepth1"] == 1
    mask_dermis_2 = adata.obs["DERdepth2"] == 1
    adata.obs.loc[mask & mask_middle_epidermis, "tissue_layer"] = "middle EPIDERMIS"
    adata.obs.loc[mask & mask_basal_epidermis, "tissue_layer"] = "basal EPIDERMIS"
    adata.obs.loc[mask & mask_upper_epidermis, "tissue_layer"] = "upper EPIDERMIS"
    adata.obs.loc[mask & mask_dermis_1, "tissue_layer"] = "DERdepth1"
    adata.obs.loc[mask & mask_dermis_2, "tissue_layer"] = "DERdepth2"
    mask = adata.obs["tissue_layer"] == "JUNCTION"
    adata.obs.loc[mask, "tissue_layer"] = 'basal EPIDERMIS'
    adata.obs["tissue_layer"] = adata.obs["tissue_layer"].astype('category')
    adata.obs["tissue_layer"] = adata.obs["tissue_layer"].cat.remove_unused_categories()

    # Rename falsely labeled spots
    if np.any((adata.obs['basal EPIDERMIS'] == 1) & (adata.obs['DERdepth1'] == 1)):
        adata.obs.loc[
            (adata.obs['basal EPIDERMIS'] == 1) & (adata.obs['DERdepth1'] == 1), 'basal EPIDERMIS'] = [0, 0]

    return adata


def rename_junction_in_spottypes(adata: anndata.AnnData) -> anndata.AnnData:
    mask = adata.obs["spot_type"] == "JUNCTION"
    mask_middle_epidermis = adata.obs["middle EPIDERMIS"] == 1
    mask_basal_epidermis = adata.obs["basal EPIDERMIS"] == 1
    mask_upper_epidermis = adata.obs["upper EPIDERMIS"] == 1
    mask_dermis = adata.obs["DERMIS"] == 1
    adata.obs.loc[mask & mask_middle_epidermis, "spot_type"] = "middle EPIDERMIS"
    adata.obs.loc[mask & mask_basal_epidermis, "spot_type"] = "basal EPIDERMIS"
    adata.obs.loc[mask & mask_upper_epidermis, "spot_type"] = "upper EPIDERMIS"
    adata.obs.loc[mask & mask_dermis, "spot_type"] = "DERMIS"
    adata.obs["spot_type"] = adata.obs["spot_type"].cat.remove_unused_categories()

    adata.obs["spot_type"] = adata.obs['spot_type'].cat.reorder_categories(
        ["upper EPIDERMIS", "middle EPIDERMIS", "basal EPIDERMIS", "DERMIS", 'MUSCLE', "VESSEL", "HAIR FOLLICLE",
         "SEBACEOUS GLAND", "SWEAT GLAND"])

    colors_spottype = dict(zip(adata.obs['spot_type'].cat.categories.to_list(), [
        "#1f77b4", "#ff7f0e", "#279e68", '#e377c2', '#8c564b', '#aa40fc', '#b5bd61', '#17becf', '#aec7e8']))

    adata.uns['spot_type_colors'] = itemgetter(*adata.obs.loc[
        mask, 'spot_type'].cat.categories.to_list())(colors_spottype)

    return adata



def add_columns_genes(adata: anndata, genes: str, label: str, count_threshold: [int, dict] = 1):
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
    if isinstance(count_threshold, int):
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
    if isinstance(count_threshold, dict):
        for disease in count_threshold.keys():
            m_disease = adata.obs['DISEASE'].str.contains(disease)
            if isinstance(genes, list):
                # Get available genes in adata
                available_genes = list(set(adata.var.index) & set(genes))

                # 2. Get bool mask for gene spot positions
                # 2.1 Find index of genes
                # varindex_genes = adata.var_names.isin(['IL32', 'DAPK1'])
                varindex_genes = np.where(np.asarray(adata.var.index)[np.newaxis, :] == np.array(
                    available_genes)[:, np.newaxis])[1]
                # 2.2 Get count matrix of genes
                if "counts" in adata.layers.keys():
                    genes_matrix = np.copy(adata.layers['counts'])[:, varindex_genes]
                else:
                    genes_matrix = np.copy(adata.X)[:, varindex_genes]
                # Identify where counts > 0
                m_genes = genes_matrix >= count_threshold[disease]
                m_array_genes = m_genes.sum(axis=1).astype(bool)

                # Add label of gene
                default_label[m_array_genes & m_disease] = label
                default_clusters[m_array_genes & m_disease] = 0
                default_counts[m_array_genes & m_disease] = np.sum(genes_matrix, axis=1)[m_array_genes & m_disease]
            else:
                # 2. Check if gene is in adata object
                if genes in adata.var_names:
                    # 2.1 Get count matrix of gene
                    if "counts" in adata.layers.keys():
                        gene_matrix = np.copy(adata.layers['counts'])[:, np.where(adata.var.index == genes)[0]]
                    else:
                        gene_matrix = np.copy(adata.X)[:, np.where(adata.var.index == genes)[0]]
                    # Identify where counts > 0
                    m_gene = gene_matrix[:, 0] >= count_threshold[disease]

                    # 2.2 Get gene spot position
                    # Add label to spots
                    default_label[m_gene & m_disease] = label
                    default_clusters[m_gene & m_disease] = 0
                    default_counts[m_gene & m_disease] = gene_matrix[:, 0][m_gene & m_disease]

    else:
        if isinstance(genes, list):
            # Get available genes in adata
            available_genes = list(set(adata.var.index) & set(genes))

            # 2. Get bool mask for gene spot positions
            # 2.1 Find index of genes
            varindex_genes = np.where(np.asarray(adata.var.index)[np.newaxis, :] == np.array(
                available_genes)[:, np.newaxis])[1]
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


def mark_spotsincluster(adata: anndata.AnnData, sub_adata: anndata.AnnData, spot_indices: [int, list],
                        obs_conditional_gene_counts: str, gene: str, radius: [int, str]):
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
