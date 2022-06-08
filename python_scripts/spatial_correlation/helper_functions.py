"""Figure parameters
    File name: helper_functions.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

import numpy as np
import pandas as pd
from collections import OrderedDict


def figure_params():
    fig_size = (6, 6)
    xy_fontsize = 16
    xy_ticks = 12
    title_fontsize = 18
    legend_fontsize = 14
    text_fontsize = 14
    fileformat = '.pdf'
    img_key = 'hires'  # 'lowres'
    return fig_size, title_fontsize, xy_fontsize, legend_fontsize, fileformat, img_key, xy_ticks, text_fontsize


def get_cropped_sampleimg(img_key):
    samples = ['P15509_1003', 'P15509_1004', 'P16357_1003', 'P16357_1019', 'P16357_1020', 'P16357_1031',
               'P16357_1032', 'P16357_1036']
    if img_key == 'lowres':
        crops_img = [(100, 400, 520, 130), (49, 349, 510, 100), (80, 380, 510, 100), (150, 450, 510, 100),
                     (150, 450, 530, 120), (135, 435, 530, 120), (150, 450, 510, 100), (160, 460, 510, 100)]
    else:
        res_proportion = 2000 / 600
        crops_img = [(100, 400, 520, 130), (49, 349, 510, 100), (80, 380, 510, 100), (150, 450, 510, 100),
                     (150, 450, 530, 120), (135, 435, 530, 120), (150, 450, 510, 100), (160, 460, 510, 100)]
        crops_img = (np.asarray(crops_img) * res_proportion).round()

    return samples, crops_img


def interface_to_epidermis(adata, tissue_layers, epidermis_layers):
    # Rename tissue region 'INTERFACE' to upper, middle or basal EPIDERMIS because some spots got both labels
    m_interface = adata.obs['JUNCTION'] == 1
    if isinstance(epidermis_layers, list):
        df_temp = adata.obs[epidermis_layers][m_interface]
    else:
        df_temp = adata.obs[[epidermis_layers]][m_interface]
    df_temp = df_temp.loc[:, df_temp.columns].replace(1, pd.Series(df_temp.columns, df_temp.columns))
    df_temp['JUNCTION'] = '0'
    for col in df_temp.columns[:-1]:
        df_temp['JUNCTION'] += df_temp[col].astype(str)
    df_temp["JUNCTION"] = df_temp.JUNCTION.str.replace('0', '')
    adata.obs[tissue_layers][m_interface] = df_temp["JUNCTION"]
    adata.obs[tissue_layers] = adata.obs[tissue_layers].cat.remove_unused_categories()

    return adata


def get_color_signaturegenes():
    # Group cytokines into Type 1-3 (Signature T-cells)
    signatures = OrderedDict()
    # publication
    signatures["IFNG"] = "#ff7f00"  # orange LICHEN
    signatures["IL13"] = "#e41a1c"  # red AE
    signatures["IL17A"] = "#377eb8"  # blue PSO
    signatures["HkG"] = '#4daf4a'  # green GAPDH
    signatures["IFNG_IL13_responder"] = "sandybrown"  # orange LICHEN
    signatures["IFNG_IL17A_responder"] = "goldenrod"  # red AE
    signatures["IL13_IFNG_responder"] = "firebrick"  # orange LICHEN
    signatures["IL13_IL17A_responder"] = "lightcoral"  # red AE
    signatures["IL17A_IFNG_responder"] = "royalblue"  # orange LICHEN
    signatures["IL17A_IL13_responder"] = "blueviolet"  # red AE

    return signatures
