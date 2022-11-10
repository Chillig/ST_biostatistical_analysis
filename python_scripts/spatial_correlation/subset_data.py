from python_scripts.spatial_correlation import tools, helper_functions as ht


def data_preparation(adata, tissue_layers, epidermis_layers, conditional_genes, conditionalgenes_responders):
    """Prepare data before applying conditional density clustering algorithm

    Parameters
    ----------
    adata : annData
    tissue_layers : list
    conditional_genes : list
    epidermis_layers : str, list
    conditionalgenes_responders : dict

    Returns
    -------

    """
    # Subset adata to tissue_types of interest: upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS'
    if tissue_layers:
        bool_col = adata.obs[tissue_layers] == 1
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