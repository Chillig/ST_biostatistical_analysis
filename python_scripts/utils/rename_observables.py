import pandas as pd


def interface_to_epidermis(adata, tissue_layers_obsname, skin_layers):
    # Rename tissue region 'INTERFACE' to upper, middle or basal EPIDERMIS because some spots got both labels
    m_interface = adata.obs['JUNCTION'] == 1
    if isinstance(skin_layers, list):
        df_temp = adata.obs[skin_layers][m_interface]
    else:
        df_temp = adata.obs[[skin_layers]][m_interface]
    df_temp = df_temp.loc[:, df_temp.columns].replace(1, pd.Series(df_temp.columns, df_temp.columns))
    df_temp['JUNCTION'] = '0'
    for col in df_temp.columns[:-1]:
        df_temp['JUNCTION'] += df_temp[col].astype(str)
    df_temp["JUNCTION"] = df_temp.JUNCTION.str.replace('0', '')
    adata.obs[tissue_layers_obsname][m_interface] = df_temp["JUNCTION"]
    adata.obs[tissue_layers_obsname] = adata.obs[tissue_layers_obsname].cat.remove_unused_categories()

    return adata
