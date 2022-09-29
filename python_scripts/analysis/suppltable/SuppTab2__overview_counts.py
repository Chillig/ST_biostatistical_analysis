import numpy as np
import pandas as pd
import os
import anndata


def create_supptable_1(unpp_st_adata: anndata, save_folder: str):
    # Suppl. Table on UMI counts and spots. This one should include the following information for each patient:
    # -          Number of spots/section
    # -          Total UMI count/section
    # -          Median UMI count/spot/section
    # -          Number IFNg+ spots
    # -          Median IFNg UMI count/IFNg+ spot
    # -          Number IL13+ spots/section
    # -          Median IL-13 UMI count/IL-13+ spot
    # -          Number IL17+ spots/section
    # -          Median IL-17 UMI count/IL-17+ spot
    suppl_table = {key: [] for key in unpp_st_adata.obs['DISEASE'].cat.categories}
    suppl_table_temp = {key: [] for key in unpp_st_adata.obs['DISEASE'].cat.categories}

    for diag in unpp_st_adata.obs['DISEASE'].cat.categories:
        unpp_st_adata_diag = unpp_st_adata[unpp_st_adata.obs['DISEASE'] == diag].copy()
        suppl_table[diag] = {key: [] for key in unpp_st_adata_diag.obs['patient'].cat.categories}

        for patient in unpp_st_adata_diag.obs['patient'].cat.categories:
            unpp_st_adata_patient = unpp_st_adata_diag[unpp_st_adata_diag.obs['patient'] == patient].copy()

            df_temp = pd.DataFrame(
                index=['Number of spots/section', 'Total UMI count/section', 'Median UMI count/spot/section',
                       'Number IFNg+ spots/section', 'Min IFNg UMI count/IFNg+ spot in section',
                       'Max IFNg UMI count/IFNg+ spot in section', 'Median IFNg UMI count/IFNg+ spot in section',
                       'Number IL13+ spots/section', 'Min IL-13 UMI count/IL-13+ spot in section',
                       'Max IL-13 UMI count/IL-13+ spot in section', 'Median IL-13 UMI count/IL-13+ spot in section',
                       'Number IL17A+ spots/section', 'Min IL-17A UMI count/IL-17A+ spot in section',
                       'Max IL-17A UMI count/IL-17A+ spot in section', 'Median IL-17A UMI count/IL-17+ spot in section'],
                columns=list(unpp_st_adata_patient.obs['specimen'].cat.categories))

            for section in unpp_st_adata_patient.obs['specimen'].cat.categories:
                unpp_st_adata_section = unpp_st_adata_patient[unpp_st_adata_patient.obs['specimen'] == section].copy()
                df_unpp_st_adata_section = unpp_st_adata_section.to_df()
                df_temp.loc['Number of spots/section', section] = unpp_st_adata_section.shape[0]
                df_temp.loc['Total UMI count/section', section] = unpp_st_adata_section.X.sum(0).sum()
                df_temp.loc['Median UMI count/spot/section', section] = df_unpp_st_adata_section.sum(1).median()
                # Get IFNG+ spots
                df_temp.loc['Number IFNg+ spots/section', section] = np.count_nonzero(df_unpp_st_adata_section['IFNG'])
                umicounts_ifng = df_unpp_st_adata_section[df_unpp_st_adata_section['IFNG'] > 0].loc[:, 'IFNG']
                df_temp.loc['Min IFNg UMI count/IFNg+ spot in section', section] = umicounts_ifng.min(axis=0)
                df_temp.loc['Max IFNg UMI count/IFNg+ spot in section', section] = umicounts_ifng.max(axis=0)
                df_temp.loc['Median IFNg UMI count/IFNg+ spot in section', section] = umicounts_ifng.median(axis=0)
                # Get IL13+ spots
                df_temp.loc['Number IL13+ spots/section', section] = np.count_nonzero(df_unpp_st_adata_section['IL13'])
                umicounts_il13 = df_unpp_st_adata_section[df_unpp_st_adata_section['IL13'] > 0].loc[:, 'IL13']
                df_temp.loc['Min IL-13 UMI count/IL-13+ spot in section', section] = umicounts_il13.min(axis=0)
                df_temp.loc['Max IL-13 UMI count/IL-13+ spot in section', section] = umicounts_il13.max(axis=0)
                df_temp.loc['Median IL-13 UMI count/IL-13+ spot in section', section] = umicounts_il13.median(axis=0)
                # Get IL17A+ spots
                df_temp.loc['Number IL17+ spots/section', section] = np.count_nonzero(df_unpp_st_adata_section['IL17A'])
                umicounts_il17a = df_unpp_st_adata_section[df_unpp_st_adata_section['IL17A'] > 0].loc[:, 'IL17A']
                df_temp.loc['Min IL-17A UMI count/IL-17A+ spot in section', section] = umicounts_il17a.min(axis=0)
                df_temp.loc['Max IL-17A UMI count/IL-17A+ spot in section', section] = umicounts_il17a.max(axis=0)
                df_temp.loc['Median IL-17A UMI count/IL-17+ spot in section', section] = umicounts_il17a.median(axis=0)

            suppl_table[diag][patient] = df_temp

        suppl_table_temp[diag] = pd.concat(suppl_table[diag], axis=1)
    pd.concat(suppl_table_temp, axis=1).to_excel(os.path.join(save_folder, 'Suppl_table_overview.xlsx'))
