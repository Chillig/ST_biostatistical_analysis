import utils


def data_preparation(adata, conditional_genes, conditionalgenes_responders):
    """Prepare data before applying conditional density clustering algorithm

    Parameters
    ----------
    adata : annData
    conditional_genes : list
    conditionalgenes_responders : dict

    Returns
    -------

    """
    # save .uns
    dict_spatial = adata.uns['spatial']

    # Subset adata to tissue_types of interest: upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS'
    adata = utils.remove_junction_tissuelayer_spottype(adata=adata)

    # Get counts of cyotkines and their responders in the EPIDERMIS
    # - distance between spots: 100 µm, spot diameter: 55 µm
    # - better: use index array
    diseases = list(adata.obs['DISEASE'].cat.categories)
    resp_label = []
    labels = ["{}_{}".format(s, "Responders") for s in conditional_genes]
    labels = labels + conditional_genes
    max_nl_counts_cytokines = dict.fromkeys(labels)
    for cyto in conditional_genes:
        adata = utils.add_columns_genes(adata=adata, genes=cyto, label=cyto, count_threshold=1)
        # # Calculate new cut-off using non lesion of same disease as reference
        # max_nl_counts_cytokines[cyto] = dict.fromkeys(diseases)
        # for disease in diseases:
        #     adata_disease = adata[adata.obs['DISEASE'] == disease].copy()
        #     max_nl_counts_cytokines[cyto][disease] = adata_disease.obs[
        #         '{}_counts'.format(cyto)][adata_disease.obs['biopsy_type'] == 'NON LESIONAL'].max()

    for ind, cyto_resps in enumerate(conditionalgenes_responders.keys()):
        resp_label.append("_".join([cyto_resps, "Responders"]))
        adata = utils.add_columns_genes(adata=adata, genes=conditionalgenes_responders[cyto_resps],
                                        label=resp_label[ind], count_threshold=1)
        # # Calculate new cut-off using non lesion of same disease as reference
        # max_nl_counts_cytokines[resp_label[ind]] = dict.fromkeys(diseases)
        # for disease in diseases:
        #     adata_disease = adata[adata.obs['DISEASE'] == disease].copy()
        #     max_nl_counts_cytokines[resp_label[ind]][disease] = adata_disease.obs[
        #         '{}_counts'.format(resp_label[ind])][adata_disease.obs['biopsy_type'] == 'NON LESIONAL'].max()
        #
        # # apply new count cut-off - overwrite old one
        # adata = tools.add_columns_genes(adata=adata, genes=conditionalgenes_responders[cyto_resps],
        #                                 label=resp_label[ind], count_threshold=max_nl_counts_cytokines[resp_label[ind]])

        # store spatial information
        adata.obs['sample'] = adata.obs['sample'].cat.remove_unused_categories()

        dict_spatial_wanted = dict(
            (k, dict_spatial[k]) for k in list(adata.obs['sample'].cat.categories) if k in dict_spatial)
        adata.uns['spatial'] = dict_spatial_wanted

    print('Cut-offs determine if a gene is expressed in comparison to non-lesion skin: ', max_nl_counts_cytokines)

    return adata
