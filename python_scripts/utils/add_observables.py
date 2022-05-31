from python_scripts.utils import helper_tools as ht

import numpy as np
from itertools import combinations


def add_lesion_metadata(adata):
    """Add metaData as observable to annData object
        -> data specific

    Parameters
    ----------
    adata : annData

    Returns
    -------

    """
    # get labels
    disease_labels, lesion_labels = ht.get_tissue_annot(adata)

    # assign biopsy type
    adata.obs['biopsy_type'] = 'Unknown'
    adata.obs['biopsy_type'] = adata.obs['biopsy_type'].astype('<U16')
    for spot_label in lesion_labels:
        m_spots = adata.obs[spot_label] == 1
        adata.obs['biopsy_type'][m_spots] = spot_label
    adata.obs['biopsy_type'] = adata.obs['biopsy_type'].astype('string').astype('category')

    return adata, lesion_labels


def add_spottypes_obs(adata):
    """Add tissue layers as observable to adata

    Parameters
    ----------
    adata : annData

    Returns
    -------

    """
    # 1. Select tissue types of interest
    # Stand: 22.03.2022
    spot_types = ['DERMIS', 'upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS', 'JUNCTION',
                  'SEBACEOUS GLAND', 'SWEAT GLAND', 'MUSCLE', 'HAIR FOLLICLE', 'VESSEL']

    adata.obs['spot_type'] = 'Unknown'
    adata.obs['spot_type'] = adata.obs['spot_type'].astype('<U64')
    for tissue in spot_types:
        m_tissue = adata.obs[tissue] == 1
        adata.obs['spot_type'][m_tissue] = tissue

    adata.obs['spot_type'] = adata.obs['spot_type'].astype('category')
    adata.obs['spot_type'] = adata.obs['spot_type'].cat.reorder_categories(
        ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS', 'JUNCTION', 'DERMIS',
         'MUSCLE', 'VESSEL', 'HAIR FOLLICLE', 'SEBACEOUS GLAND', 'SWEAT GLAND', 'Unknown'])

    return adata


def add_tissuelayers_obs(adata):
    """Add tissue layers as observable to adata

    Parameters
    ----------
    adata : annData

    Returns
    -------

    """
    # 1. Select tissue types of interest
    # Stand: 22.03.2022
    tissue_types = ['DERdepth1', 'DERdepth2', 'DERdepth3', 'DERdepth4', 'DERdepth5', 'DERdepth6', 'DERdepth7',
                    'upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS', 'JUNCTION']

    adata.obs['tissue_layer'] = 'Unknown'
    adata.obs['tissue_layer'] = adata.obs['tissue_layer'].astype('<U64')
    for tissue in tissue_types:
        m_tissue = adata.obs[tissue] == 1
        adata.obs['tissue_layer'][m_tissue] = tissue

    adata.obs['tissue_layer'] = adata.obs['tissue_layer'].astype('category')
    adata.obs['tissue_layer'] = adata.obs['tissue_layer'].cat.reorder_categories(
        ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS', 'JUNCTION', 'DERdepth1', 'DERdepth2', 'DERdepth3',
         'DERdepth4', 'DERdepth5', 'DERdepth6', 'DERdepth7', 'Unknown'])

    return adata


def add_disease_healthy_obs(adata):
    adata.obs['healthy_disease'] = 'NON LESIONAL'
    adata.obs['healthy_disease'] = adata.obs['healthy_disease'].astype('<U32')
    m_disease = adata.obs['NON LESIONAL'] == 0
    adata.obs['healthy_disease'][m_disease] = adata.obs['DISEASE'][m_disease]

    return adata


def convert_variable_to_observable(adata, gene_names, task='cell_gene', label='celltype', condition=None):
    """Add observable column to annData by gene(s) of interest

    Parameters
    ----------
    adata : annData
    gene_names : dict or list or str
        gene(s) of interest
    task : str
        which task should be done
    label : str
        which observable name to add to annData.obs
    condition : 'list' ['numpy']
        contains functions

    Returns
    -------

    """
    if task == 'cell_gene':
        if isinstance(gene_names, list):
            obs_name = []
            for gen in gene_names:
                adata, observable = _mark_cell_by_gene(adata=adata, gene=gen)
                obs_name.append(observable)
        else:
            adata, obs_name = _mark_cell_by_gene(adata=adata, gene=gene_names)
    elif task == 'dge_prep':
        adata, obs_name = _get_conditions_dge(adata=adata, genes=gene_names)
    elif task == 'annotate_cells':
        adata, obs_name = _mark_cells_by_markergene(adata=adata, genes_dict=gene_names, label=label,
                                                    condition=condition)
    else:
        obs_name = ""

    return adata, obs_name


def _mark_cell_by_gene(adata, gene):
    """
    Mark all cells having at least 1 UMI-count of gene as gene+ cells and all others as gene- cells
    1. Adds gene+ and gene- obs. column
    2. Adds counts for specific gene in gene+ cell and counts of gene- cells to adata.obs

    :param adata: [annData]
    :param gene: [string] gene name
    :return:
    """
    # Get column index with gene name
    # Get rows where gene counts > 0
    # Sum up counts
    # Get corresponding barcodes
    # get corresponding sample id
    # get corresponding index

    if gene in adata.var.index:
        # get counts
        if "counts" in adata.layers.keys():
            counts_cyto = np.copy(adata.layers["counts"])[:, np.where(adata.var.index == gene)[0]]
            summed_counts_bc = np.copy(adata.layers["counts"]).sum(axis=1)
        else:
            counts_cyto = np.copy(adata.X)[:, np.where(adata.var.index == gene)[0]]
            summed_counts_bc = np.copy(adata.X).sum(axis=1)
        counts_cyto_gene = counts_cyto[:, 0]
        # create mask
        m_cyto = counts_cyto_gene > 0

        # Add label column to adata.obs
        obs_name = "_".join(["cytokine", gene])
        adata.obs[obs_name] = gene
        adata.obs[obs_name].loc[~m_cyto] = 'Others'
        adata.obs[obs_name] = adata.obs[obs_name].astype('category')

        # add counts column to adata.obs which contains counts of cytokine and others
        adata.obs["_".join([gene, "counts"])] = counts_cyto_gene
        # find all rows in matrix which dont belong to the cytokine
        adata.obs["_".join([gene, "counts"])].loc[~m_cyto] = summed_counts_bc[~m_cyto]
    else:
        obs_name = ""

    return adata, obs_name


def _get_conditions_dge(adata, genes):
    observables = []
    # Get available genes
    available_genes = list(set(adata.var.index) & set(genes))

    # Get index of cytokines and other genes; the position of the index = position of cytokine in available_cytos list
    varindex_cyto_genes = np.where(adata.var.index[np.newaxis, :] == np.array(available_genes)[:, np.newaxis])[1]
    # get counts
    if "counts" in adata.layers.keys():
        counts_cyto = np.copy(adata.layers["counts"])[:, varindex_cyto_genes]
        summed_counts_bc = np.sum(adata.layers['counts'], axis=1)
    else:
        counts_cyto = np.copy(adata.X)[:, varindex_cyto_genes]
        summed_counts_bc = np.sum(adata.X, axis=1)
    # create mask
    m_cyto = counts_cyto > 0

    # Loop through each cytokine
    m_cytokines = np.ones(m_cyto[:, 0].shape[0]).astype(bool)
    for ind, cyto in enumerate(available_genes):
        # get combined mask using the or operator
        m_cytokines = m_cytokines + m_cyto[:, ind]

        # Add label column to adata.obs
        obs_name = "_".join([cyto, 'others'])
        adata.obs[obs_name] = cyto
        adata.obs[obs_name].loc[~m_cyto[:, ind]] = 'Others'

        obs_name = "_".join([cyto, 'cyto_posneg'])
        adata.obs[obs_name] = "condition 1"
        adata.obs[obs_name].loc[~m_cyto[:, ind]] = "condition 2"

        # add counts column to adata.obs which contains counts of cytokine and others
        obs_name = "_".join(["counts_1", cyto])
        adata.obs[obs_name] = 1
        adata.obs[obs_name].loc[m_cyto[:, ind]] = counts_cyto[m_cyto[:, ind]].T[ind]
        # find all rows in matrix which dont belong to the cytokine
        adata.obs[obs_name].loc[~m_cyto[:, ind]] = summed_counts_bc[~m_cyto[:, ind]]

        observables.append(obs_name)

    return adata, observables


def _mark_cells_by_markergene(adata, genes_dict, label, condition):
    """
    Attention: This function already assumes that all genes are in adata object!!
    Mark cells or spots by a marker gene.
    If a spot contains more than one marker gene it is marked as double, triple ... positive cell or spot

    :param adata: [annData]
    :param genes_dict: [dict]
    :param label: [string] label or observable name
    :param condition: [numpy.operation]
        for AND operation on a matrix use np.all for OR operation use np.sum or np.any
    :return:
    """
    # Get CD4 and CD8 cells which contain CD4+ or CD8A+/CD8B+

    obs_counts = "_".join([label, 'counts'])
    obs_label = "_".join([label, 'others'])

    # mark cells which don't contain those genes as others
    default_label = np.array(['Others'] * adata.shape[0], dtype='<U32')

    # initialize mask
    mask_matrix_genes = np.zeros(shape=(adata.shape[0], len(genes_dict.keys())))
    index_genes = []

    # get counts
    if "counts" in adata.layers.keys():
        default_counts = np.copy(adata.layers['counts'].sum(axis=1))
    else:
        default_counts = np.copy(adata.X.sum(axis=1))

    for ind, cell in enumerate(genes_dict.keys()):
        # First check if genes are in data set
        if isinstance(genes_dict[cell], str):
            available_genes = list(set(adata.var.index) & {genes_dict[cell]})
        elif len(genes_dict[cell]) == 0:
            print('Not yet implemented..')  # TODO case 0 gene .. should not happen as you as a user provide the names..
        else:
            available_genes = list(set(adata.var.index) & set(genes_dict[cell]))

        if len(available_genes) > 0:
            # Get index of genes and other genes; the position of the index = position of genes in available_genes list
            varindex_genes = np.where(adata.var.index[np.newaxis, :] == np.array(available_genes)[:, np.newaxis])[1]
            index_genes.extend(varindex_genes)
            # get counts
            if "counts" in adata.layers.keys():
                counts_genes = adata.layers["counts"][:, varindex_genes].copy()
            else:
                counts_genes = adata.X[:, varindex_genes].copy()

            # create mask
            if counts_genes.shape[1] > 1:
                m_genes = counts_genes > 0
                # Identify where counts > 0:
                #   for AND operation on a matrix use np.all (bool array) for
                #   OR operation use np.sum (integer array) or np.any (bool array)
                m_array_genes = condition[ind](m_genes, axis=1)
                # get counts
                default_counts[m_array_genes] = np.sum(counts_genes, axis=1)[m_array_genes]
            else:
                m_array_genes = counts_genes[:, 0] > 0
                counts_genes = counts_genes[:, 0]
                default_counts[m_array_genes] = counts_genes[m_array_genes]

            # # Add label of cell type
            default_label[m_array_genes] = cell
            mask_matrix_genes[:, ind] = m_array_genes

    # apply double positive label to those and get counts
    dp_genes = np.asarray(adata.var_names[np.unique(index_genes)])
    # get position of double positive cells
    if len(dp_genes) == 2:
        indices = np.argwhere(np.all(mask_matrix_genes == 1, axis=1))[:, 0]
        default_label[indices] = " & ".join(dp_genes)
    else:
        # Get all possible combinations for length 2 to number of double positive (dp) genes
        comb_dp_genes = []
        for n in range(2, len(dp_genes) + 1):
            comb_dp_genes.extend([i for i in combinations(dp_genes, n)])
        indices = []
        # keys as array
        dictkey_array = np.asarray(list(genes_dict.keys()))
        for combi_genes in comb_dp_genes:
            # find out the position of the genes in the combination in the dictionary
            index_combgenes = np.where(dictkey_array[np.newaxis, :] == np.array(list(combi_genes))[:, np.newaxis])[1]
            # read out indices and save joint label
            temp_indices = np.argwhere(np.all(mask_matrix_genes[:, index_combgenes] == 1, axis=1))[:, 0]
            default_label[temp_indices] = " & ".join(list(combi_genes))
            indices.extend(temp_indices)

    # get counts
    if "counts" in adata.layers.keys():
        counts_genes = np.copy(adata.layers["counts"])[:, index_genes]
    else:
        counts_genes = np.copy(adata.X)[:, index_genes]
    default_counts[indices] = np.sum(counts_genes[indices], axis=1)

    # 4. Add observables to adata
    adata.obs[obs_label] = default_label
    adata.obs[obs_counts] = default_counts

    return adata, obs_label


def add_columns_genes(adata, genes, label):
    """
    Add columns of counts, label, and clusters to adata object beginning with term label
    gene(s) label = 0
    others label = 1

    :param adata: [annData] object which contains all information about the experiment
    :param genes: [string or list] it is possible to only provide a gene name or a list of gene names
    :param label: [string] label for new column
    :return: [annData]
    """

    # name columns
    obs_counts = "_".join([label, 'counts'])
    obs_label = "_".join([label, 'label'])
    obs_clusters = "_".join([label, 'clusters'])

    # set default counts to 0
    default_counts = np.array([0] * adata.shape[0])
    # set default label "Others"
    default_label = np.array(['Others'] * adata.shape[0], dtype='<U24')
    default_clusters = np.array([1] * adata.shape[0])

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
        m_genes = genes_matrix > 0
        m_array_genes = m_genes.sum(axis=1).astype(bool)

        # Add label of responder
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
            m_gene = gene_matrix[:, 0] > 0

            # 2.2 Get gene spot position
            # Add label cytokine positive spots
            default_label[m_gene] = label
            default_clusters[m_gene] = 0
            default_counts[m_gene] = gene_matrix[:, 0][m_gene]

    # 3. Add observables to adata
    adata.obs[obs_clusters] = default_clusters
    adata.obs[obs_label] = default_label
    adata.obs[obs_counts] = default_counts

    # 4. Convert to categorical
    adata.obs[obs_clusters] = adata.obs[obs_clusters].astype('category')
    adata.obs[obs_label] = adata.obs[obs_label].astype('category')

    return adata
