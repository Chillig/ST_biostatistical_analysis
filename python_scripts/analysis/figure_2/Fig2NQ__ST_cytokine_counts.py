import pandas as pd
import os


def get_per_diagnosis_cytokine_counts(unpp_st_adata, save_folder, biopsy_type):
    # Read out only lesion or Non-lesion biopsies
    unpp_st_adata = unpp_st_adata[unpp_st_adata.obs['biopsy_type'] == biopsy_type].copy()

    cytokines = ['IL17A', 'IL13', 'IFNG', 'TNF', 'IL17F', 'IL21', 'IL22', 'IL10', 'IL4']
    dict_cytokine_counts = {}
    for cyto in cytokines:
        unpp_st_adata_temp = unpp_st_adata[:, unpp_st_adata.var_names == cyto].copy()

        df_cytokine = pd.DataFrame()
        for diag in ['Pso', 'AD', 'LP']:
            unpp_st_adata_temp_diag = unpp_st_adata_temp[unpp_st_adata_temp.obs['DISEASE'] == diag].copy()

            cyto_counts = []
            for specimen in unpp_st_adata_temp_diag.obs['specimen'].cat.categories:
                cyto_counts.append(unpp_st_adata_temp_diag[
                                       unpp_st_adata_temp_diag.obs['specimen'] == specimen].to_df().sum().values[0])
            df_temp = pd.DataFrame({diag: cyto_counts})
            df_cytokine = pd.concat([df_cytokine, df_temp], ignore_index=True, axis=1)
        df_cytokine = df_cytokine.rename({0: 'Pso', 1: 'AD', 2: 'LP'}, axis='columns')

        dict_cytokine_counts[cyto] = df_cytokine
    pd.concat(dict_cytokine_counts, axis=1).to_excel(
        os.path.join(save_folder, '{}_Tcell_cytokinecounts.xlsx'.format(biopsy_type)))


def get_cytokinecounts_per_diagnosis(unpp_st_adata, save_folder, biopsy_type):
    # Read out only lesion or Non-lesion biopsies
    unpp_st_adata = unpp_st_adata[unpp_st_adata.obs['biopsy_type'] == biopsy_type].copy()

    cytokines = ['IL17A', 'IL13', 'IFNG', 'TNF', 'IL17F', 'IL21', 'IL22', 'IL10', 'IL4']
    df_cytopercentage = pd.DataFrame(index=cytokines)
    for diag in ['Pso', 'AD', 'LP']:
        unpp_st_adata_temp = unpp_st_adata[unpp_st_adata.obs['DISEASE'] == diag].copy()
        df_cytocounts = pd.DataFrame(index=cytokines)

        for specimen in unpp_st_adata_temp.obs['specimen'].cat.categories:
            cyto_counts = unpp_st_adata_temp[
                unpp_st_adata_temp.obs['specimen'] == specimen].to_df()[cytokines].sum(axis=0)
            cyto_counts = pd.DataFrame({specimen: cyto_counts.values}, index=cyto_counts.index)
            # if cyto_counts[specimen].sum() != 0:
            df_cytocounts = pd.concat([df_cytocounts, cyto_counts], axis=1)

        df_cytocounts.to_excel(os.path.join(save_folder, '{}_{}_cytokinecounts.xlsx'.format(biopsy_type, diag)))

        # Fig2Q
        df_diagcyto_summed = df_cytocounts.transpose().sum(axis=0)
        # total_counts_disease = df_diagcyto_summed.sum()
        # df_temp = df_diagcyto_summed / total_counts_disease
        df_temp = pd.DataFrame({diag: df_diagcyto_summed.values}, index=df_diagcyto_summed.index)

        df_cytopercentage = pd.concat([df_cytopercentage, df_temp], axis=1)

    df_cytopercentage.to_excel(os.path.join(save_folder, '{}_Disease_cytokinecounts.xlsx'.format(biopsy_type)))
