"""Figure parameters
    File name: helper_functions.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

import numpy as np


def figure_params():
    fig_size = (8, 8)
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
