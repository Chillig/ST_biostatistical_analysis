import pandas as pd


def interface_to_epidermis(adata, tissue_layers):
    # Rename tissue region 'INTERFACE' to upper, middle or basal EPIDERMIS because some spots got both labels
    m_interface = adata.obs['tissue_type'] == 'INTERFACE'
    if isinstance(tissue_layers, list):
        df_temp = adata.obs[tissue_layers][m_interface]
    else:
        df_temp = adata.obs[[tissue_layers]][m_interface]
    df_temp = df_temp.loc[:, df_temp.columns].replace(1, pd.Series(df_temp.columns, df_temp.columns))
    df_temp['interface'] = '0'
    for col in df_temp.columns[:-1]:
        df_temp['interface'] += df_temp[col].astype(str)
    df_temp["interface"] = df_temp.interface.str.replace('0', '')
    adata.obs['tissue_type'][m_interface] = df_temp["interface"]
    adata.obs['tissue_type'] = adata.obs['tissue_type'].cat.remove_unused_categories()

    return adata
