from python_scripts.pre_processing.plot_functions import helper_functions as hf

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os

fig_size, title_fontsize, subtitle_fontsize, xy_fontsize, xy_ticks, legend_fontsize, text_fontsize, fileformat = \
    hf.plot_properties()


def plot_greyspots_image(configs, adata, save_folder, label):
    """Plot spots in grey on top of H&E slide

    Parameters
    ----------
    configs : configparser
    adata : annData
    save_folder : str
    label : None, str

    Returns
    -------

    """
    # Size of count spots
    size = 0.9
    samples = np.unique(adata.obs['sample'].values)
    for ind, sample in enumerate(samples):
        subadata = adata[adata.obs['sample'] == sample].copy()
        fig, ax = plt.subplots(figsize=fig_size)
        sc.pl.spatial(subadata, color=label, size=size, library_id=sample, ax=ax,
                      show=False, img_key=configs['visualisation_options']['image_res'])

        # Invert both axis due to flipped and mirrored images
        ax.invert_xaxis()
        ax.invert_yaxis()

        # Legend: outside of axis 1.45
        leg = ax.legend(bbox_to_anchor=(1.45, 0.6), ncol=1, fontsize=legend_fontsize)
        leg.get_frame().set_linewidth(0.0)

        plt.tight_layout()
        if configs.getboolean('preprocessing', 'read_raw_matrix'):
            plt.savefig(os.path.join(save_folder, "_".join(["raw_H&E_image", sample, fileformat])),
                        bbox_inches='tight', bbox_extra_artists=(leg,))
        else:
            plt.savefig(os.path.join(save_folder, "_".join(["H&E_image", sample, fileformat])),
                        bbox_inches='tight', bbox_extra_artists=(leg,))
        plt.close()
