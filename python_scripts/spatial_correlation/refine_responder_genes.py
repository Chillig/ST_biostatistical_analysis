import scanpy as sc
import os
import pandas as pd


def rank_densitycluster_vs_others_genes(adata, cytokine_responders, save_folder, radius):
    for cyto in cytokine_responders.keys():
        adata.obs['{}_in_sdcc'.format(cyto)] = adata.obs['{}_in_sdcc'.format(cyto)].astype(str)
        adata.obs['{}_in_sdcc_merged'.format(cyto)] = adata.obs['{}_in_sdcc'.format(cyto)].copy()
        # decide if in cluster or no
        adata.obs['{}_in_sdcc_merged'.format(cyto)][adata.obs['{}_in_sdcc_merged'.format(cyto)] == 2] = 1
        adata.obs['{}_in_sdcc_merged'.format(cyto)] = adata.obs[
            '{}_in_sdcc_merged'.format(cyto)].astype(str).astype('category')
        sc.tl.rank_genes_groups(adata, '{}_in_sdcc_merged'.format(cyto), method='wilcoxon',
                                key_added="{}_in_sdcc_merged_wilcoxon".format(cyto),
                                use_raw=True, n_genes=adata.shape[1], corr_method='benjamini-hochberg')
        # save as pandas dataframe with p-values
        cluster_algo = '{}_in_sdcc_merged'.format(cyto)
        ranked_genes_list = adata.uns[cluster_algo + '_wilcoxon']
        name_csv_file = os.path.join(save_folder, cluster_algo + '_clustered_gene_list.csv')
        groups = ranked_genes_list['names'].dtype.names
        gene_p_value_df = pd.DataFrame({group + '_' + key: ranked_genes_list[key][group]
                                        for group in groups
                                        for key in ['names', 'scores', 'logfoldchanges',
                                                    'pvals', 'pvals_adj']})

        # Add gene_name (ENSG) to dataframe
        ind = 0
        first_cluster = min(gene_p_value_df.columns.str.split('_', 1).str[0].astype(int))
        for cluster in range(first_cluster, len(adata.obs[cluster_algo].cat.categories) + first_cluster):
            gene_p_value_df.insert(ind, '{}_ensemblID'.format(cluster),
                                   adata.var.loc[gene_p_value_df['{}_names'.format(cluster)], 'gene_id'].values)
            ind += 6  # 6: because we have 5 columns already per cluster score, names, ..

        # Plot Ranking
        # sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key="IL17A_in_sdcc_merged_wilcoxon")
        # sc.pl.rank_genes_groups_heatmap(adata, groups=['0', '1'], n_genes=5, key='IL17A_in_sdcc_merged_wilcoxon',
        #                                 show_gene_labels=True, show=False)
        # Volcano plot -> R
        # Save df
        gene_p_value_df.to_csv(os.path.join(save_folder, '{}_in_sdcc_merged_wilcoxon__radius{}.csv'.format(
            cyto, radius)))


def rank_densitycluster_genes(adata, cytokine_responders, save_folder, radius):
    # Compare nn spots against cyto+ spots
    for cyto in cytokine_responders.keys():
        adata.obs['{}_in_sdcc'.format(cyto)] = adata.obs['{}_in_sdcc'.format(cyto)].astype(str).astype('category')
        # Read out spots belonging to a cluster -> label 1 or 2
        adata_cyto = adata[adata.obs['{}_in_sdcc'.format(cyto)] != 0].copy()
        adata_cyto.obs['{}_in_sdcc'.format(cyto)] = adata_cyto.obs['{}_in_sdcc'.format(cyto)].astype(str)
        # Workaround for scanpy: rename 1 to 0 and 2 to 1
        adata_cyto.obs.loc[adata_cyto.obs['{}_in_sdcc'.format(cyto)] == '1', '{}_in_sdcc'.format(cyto)] = '0'
        adata_cyto.obs.loc[adata_cyto.obs['{}_in_sdcc'.format(cyto)] == '2', '{}_in_sdcc'.format(cyto)] = '1'
        sc.tl.rank_genes_groups(adata_cyto, '{}_in_sdcc'.format(cyto), method='wilcoxon',
                                key_added="{}_in_sdcc_wilcoxon".format(cyto),
                                use_raw=True, n_genes=adata.shape[1], corr_method='benjamini-hochberg')
        # save as pandas dataframe with p-values
        cluster_algo = '{}_in_sdcc'.format(cyto)
        ranked_genes_list = adata_cyto.uns[cluster_algo + '_wilcoxon']
        groups = ranked_genes_list['names'].dtype.names
        gene_p_value_df = pd.DataFrame({group + '_' + key: ranked_genes_list[key][group]
                                        for group in groups
                                        for key in ['names', 'scores', 'logfoldchanges',
                                                    'pvals', 'pvals_adj']})

        # Add gene_name (ENSG) to dataframe
        ind = 0
        first_cluster = min(gene_p_value_df.columns.str.split('_', 1).str[0].astype(int))
        for cluster in range(first_cluster, len(adata_cyto.obs[cluster_algo].cat.categories) + first_cluster):
            gene_p_value_df.insert(ind, '{}_ensemblID'.format(cluster),
                                   adata_cyto.var.loc[gene_p_value_df['{}_names'.format(cluster)], 'gene_id'].values)
            ind += 6  # 6: because we have 5 columns already per cluster score, names, ..

        # Plot Ranking
        # sc.pl.rank_genes_groups(adata_cyto, n_genes=25, sharey=False, key="IL17A_in_sdcc_wilcoxon")
        # sc.pl.rank_genes_groups_heatmap(adata, groups=['0', '1'], n_genes=5, key='IL17A_in_sdcc_merged_wilcoxon',
        #                                 show_gene_labels=True, show=False)
        # Volcano plot -> R
        # Save df
        gene_p_value_df.to_csv(os.path.join(save_folder, '{}_in_sdcc_wilcoxon__radius{}.csv'.format(cyto, radius)))


def rank_cyto_vs_others_genes(adata, cytokine_responders, save_folder, radius):
    for cyto in cytokine_responders.keys():
        adata.obs['{}_in_sdcc_r{}'.format(cyto, radius)] = adata.obs[
            '{}_in_sdcc_r{}'.format(cyto, radius)].astype(str).astype('category')
        # remove nn spots
        adata_wonn = adata[~adata.obs['{}_in_sdcc_r{}'.format(cyto, radius)].str.contains('2')].copy()
        sc.tl.rank_genes_groups(adata_wonn, '{}_in_sdcc_r{}'.format(cyto, radius), method='wilcoxon',
                                key_added="{}_in_sdcc_r{}_wilcoxon".format(cyto, radius),
                                use_raw=True, n_genes=adata_wonn.shape[1], corr_method='benjamini-hochberg')
        # save as pandas dataframe with p-values
        cluster_algo = '{}_in_sdcc_r{}'.format(cyto, radius)
        ranked_genes_list = adata_wonn.uns[cluster_algo + '_wilcoxon']
        name_csv_file = os.path.join(save_folder, cluster_algo + '_clustered_gene_list.csv')
        groups = ranked_genes_list['names'].dtype.names
        gene_p_value_df = pd.DataFrame({group + '_' + key: ranked_genes_list[key][group]
                                        for group in groups
                                        for key in ['names', 'scores', 'logfoldchanges',
                                                    'pvals', 'pvals_adj']})

        # Add gene_name (ENSG) to dataframe
        ind = 0
        first_cluster = min(gene_p_value_df.columns.str.split('_', 1).str[0].astype(int))
        for cluster in range(first_cluster, len(adata_wonn.obs[cluster_algo].cat.categories) + first_cluster):
            gene_p_value_df.insert(ind, '{}_ensemblID'.format(cluster),
                                   adata_wonn.var.loc[gene_p_value_df['{}_names'.format(cluster)], 'gene_id'].values)
            ind += 6  # 6: because we have 5 columns already per cluster score, names, ..

        # Plot Ranking
        # sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key="IL17A_in_sdcc_merged_wilcoxon")
        # sc.pl.rank_genes_groups_heatmap(adata, groups=['0', '1'], n_genes=5, key='IL17A_in_sdcc_merged_wilcoxon',
        #                                 show_gene_labels=True, show=False)
        # Volcano plot -> R
        # Save df
        gene_p_value_df.to_csv(os.path.join(save_folder, '{}_in_sdcc_wilcoxon__radius{}.csv'.format(cyto, radius)))
        # sc.write(os.path.join(save_folder, '{}_in_sdcc__radius{}.h5'.format(cyto, radius)), adata)
