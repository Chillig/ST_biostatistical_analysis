"""Prepare, subset, and add metaData to annData object
    File name: utils.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

import numpy as np
import copy
import math
import pandas as pd
from operator import itemgetter


def prepare_rawadata(adata):
    """Add metadata to anndata object: n_counts, log_counts, n_genes, mt_frac, patient, tissue_type

    Parameters
    ----------
    adata : annData

    Returns
    -------

    """
    adata.obs['n_counts'] = adata.X.sum(1)
    # number of counts per spot in log
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    # number of genes per spot
    adata.obs['n_genes'] = (adata.X > 0).sum(1)
    # Check if all genes starting with MT- are actually Mitochondrial genes..
    # MT genes: ['MT-ND1', 'MT-ND2', 'MT-CO1', 'MT-CO2', 'MT-ATP8', 'MT-ATP6', 'MT-CO3',
    # 'MT-ND3', 'MT-ND4L', 'MT-ND4', 'MT-ND5', 'MT-ND6', 'MT-CYB']
    # MT-ND: mitochondrially encoded NADH dehydrogenase
    # MT-CO: mitochondrially encoded cytochrome c oxidase
    # MT-ATP: mitochondrially encoded ATP synthase
    mt_gene_mask = [gene.startswith('MT-') for gene in adata.var_names]
    adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1) / adata.obs['n_counts']

    # 1. Add meta data like which samples belong to which donor (optional)
    if "patient" not in adata.obs_keys():
        adata, tissue_cell_labels, disease_labels, lesion_labels = add_metadata(adata)
        # 1.2 Remove spots having no tissue/cell labels (since 06.10.2020)
        adata = adata[np.where(adata.obs[tissue_cell_labels].to_numpy().any(axis=1))[0]]

    adata = add_tissue_obs(adata)

    return adata


def add_metadata(adata):
    """Assign batch and patient as co-viariate to each spot (-> important for e.g. batch correction, DGE analysis ..)

    Parameters
    ----------
    adata : annData

    Returns
    -------

    """

    # get batches (object slide = batch)
    df_batches = map_sample_batch_list(adata=adata,
                                       num_samples_patient=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4])
    # assign batch to each spot
    batches = []
    for project in df_batches:
        for ind, s_c in enumerate(df_batches[project]['sample']):
            no_spots = adata[adata.obs['sample'] == s_c].shape[0]
            batches.extend(np.ones(no_spots) * df_batches[project]['batch'][ind])
    # need information about how often samples were found
    adata.obs['batch'] = batches
    adata.obs['batch'] = adata.obs['batch'].astype('int').astype('category')

    # get batches assigned to each patient (currently we have 2 to 4 samples per patient)
    # last four biopsies of last donor are from two different time points
    df_patients = map_sample_batch_list(adata=adata,
                                        num_samples_patient=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 4])
    # assign donor to each spot
    donors = []
    for project in df_patients:
        for ind, s_c in enumerate(df_patients[project]['sample']):
            no_spots = adata[adata.obs['sample'] == s_c].shape[0]
            donors.extend(np.ones(no_spots) * df_patients[project]['batch'][ind])
    # need information about how often samples were found
    adata.obs['patient'] = donors
    adata.obs['patient'] = adata.obs['patient'].astype('int').astype('category')

    # assign Diagnostics
    tissue_cell_labels, disease_labels, lesion_labels = get_tissue_annot(adata)
    spot_disease = []
    for spot_label in disease_labels:
        spots = adata[adata.obs[spot_label] == 1]
        spot_disease.extend([spot_label] * spots.shape[0])
    adata.obs['disease'] = spot_disease
    adata.obs['disease'] = adata.obs['disease'].astype('string').astype('category')

    # assign Diagnostics
    spot_biopsy = []
    for spot_label in lesion_labels:
        spots = adata[adata.obs[spot_label] == 1]
        spot_biopsy.extend([spot_label] * spots.shape[0])
    adata.obs['biopsy_type'] = spot_biopsy
    adata.obs['biopsy_type'] = adata.obs['biopsy_type'].astype('string').astype('category')

    return adata, tissue_cell_labels, disease_labels, lesion_labels


def convert_categories_cytokines_responders_others(adata, cyto_responder_genes):
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


def add_columns_genes(adata, genes, label, count_threshold=1):
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


def add_tissue_obs(adata):
    """Add tissue layers as observable to adata

    Parameters
    ----------
    adata : annData

    Returns
    -------

    """
    # 1 Select tissue types of interest
    tissue_types = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS',
                    'DERdepth1', 'DERdepth2', 'DERdepth3', 'DERdepth4', 'DERdepth5', 'DERdepth6', 'DERdepth7',
                    'INTERFACE']

    adata.obs['tissue_type'] = 'Unknown'
    adata.obs['tissue_type'] = adata.obs['tissue_type'].astype('<U64')
    for tissue in tissue_types:
        m_tissue = adata.obs[tissue] == 1
        adata.obs['tissue_type'][m_tissue] = tissue

    adata.obs['tissue_type'] = adata.obs['tissue_type'].astype('category')
    return adata


def map_sample_batch_list(adata, num_samples_patient):
    """

    Parameters
    ----------
    adata : annData
    num_samples_patient : list
        [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2]

    Returns
    -------

    """

    def order_of_magnitude(number):
        return math.floor(math.log(number, 10))

    metadata_batch = pd.DataFrame(columns=["sample", "batch"])
    metadata = dict()
    list_project_samples = dict()
    for ind, project in enumerate(np.unique(adata.obs['project'])):
        list_project_samples[project] = np.unique(adata.obs['sample'][adata.obs['project'] == project])
        sample_names = [int(i.split('_')[1]) for i in list_project_samples[project]]

        temp_num_oom = 10 ** order_of_magnitude(np.amax(sample_names))
        max_samplenumber = np.amax(sample_names) - temp_num_oom

        rest = max_samplenumber % num_samples_patient[ind]
        if rest != 0:
            max_samplenumber = max_samplenumber + rest

        # determine no batches
        num_batches = 0
        list_batches = num_samples_patient.copy()
        number_samples = copy.copy(max_samplenumber)
        while True:
            list_batches.remove(num_samples_patient[num_batches])
            if ((number_samples - num_samples_patient[num_batches]) == 0) | \
                    (num_batches == (len(num_samples_patient) - 1)):
                num_batches += 1
                break
            number_samples = number_samples - num_samples_patient[num_batches]
            num_batches += 1

        batches = []
        for i_c in range(num_batches):
            batches.extend([i_c + 1 + ind] * num_samples_patient[i_c])

        array_sampleids = np.arange(1 + temp_num_oom, temp_num_oom + max_samplenumber + 1)

        for i_c in range(len(array_sampleids)):
            metadata_batch.loc[i_c] = ["_".join([project, str(array_sampleids[i_c])]), batches[i_c]]

        metadata[project] = pd.DataFrame(columns=["sample", "batch"])
        for i_c in range(len(sample_names)):
            sample_batch = metadata_batch[(metadata_batch ==
                                           list_project_samples[project][i_c]).any(1)].stack().unique()
            metadata[project].loc[i_c] = [sample_batch[0], sample_batch[1]]
        num_samples_patient = list_batches

    return metadata


def get_tissue_annot(adata):
    """Get tissue annotations saved in adata object

    Parameters
    ----------
    adata

    Returns
    -------

    """
    # Use annotations from pathologist instead of clusters
    obs_keys = list(adata.obs_keys())
    # get all manual annotations by extracting all keys with upper case characters
    annotations = [char for char in obs_keys if any(c.isupper() for c in char)]

    if "ANNOTATOR" in annotations:
        annotations.remove("ANNOTATOR")

    # split annotations by disease and tissue / cell types
    elements = ["LESONAL", "NON LESIONAL"]
    try:
        index = []
        for el in elements:
            index.append(annotations.index(el))
        target_index = np.amax(index) + 1
        disease_index = np.amin(index)
    except ValueError:
        target_index = None
        print("REMINDER: Insert 'LESONAL' and 'NON LESIONAL' in your excel file")

    tissue_cell_labels = annotations[target_index:]
    disease_labels = np.array(annotations[:disease_index])
    lesion_labels = np.array(itemgetter(*index)(annotations))

    # remove cell cycle annotations from tissue_cell_labels list
    for cc in ["G1_score", "G2M_score", "S_score", "M_score"]:
        try:
            tissue_cell_labels.remove(cc)
        except ValueError:
            continue

    return tissue_cell_labels, disease_labels, lesion_labels
