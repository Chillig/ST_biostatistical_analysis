import numpy as np


def set_title(title=None, sup_title=None, axes=None, figure=None):
    if not sup_title:
        axes.set_title(title, fontsize=_title_fontsize)
    elif not title:
        figure.suptitle(sup_title, fontsize=_subtitle_fontsize)
    elif sup_title and title:
        axes.set_title(title, fontsize=_title_fontsize)
        figure.suptitle(sup_title, fontsize=_subtitle_fontsize)


def set_axes_label(x_label, y_label, axes):
    axes.set_xlabel(x_label)
    axes.set_ylabel(y_label)
    axes.yaxis.label.set_size(_xy_fontsize)
    axes.xaxis.label.set_size(_xy_fontsize)


def calc_sigma_containment_radius(adata, sigma_value):
    sigma_1_5_containment_radius = np.percentile(adata.obs['n_genes'], sigma_value)
    return sigma_1_5_containment_radius


def place_cut_txtbox(threshold_value, ax, y_pos_txtbox=500):
    # calculate sigma containment radius
    ax.axvline(x=threshold_value, linewidth=1, color='k', linestyle='--', zorder=3)
    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
    ax.text(threshold_value, y_pos_txtbox, str(threshold_value), ha="center", va="center", size=10, bbox=bbox_props)


def plot_properties():
    # figure properties
    fig_size = (8, 8)
    xy_fontsize = 16
    xy_ticks = 12
    title_fontsize = 18
    subtitle_fontsize = 16
    legend_fontsize = 12
    text_fontsize = 12
    fileformat = '.pdf'

    return \
        fig_size, title_fontsize, subtitle_fontsize, xy_fontsize, xy_ticks, legend_fontsize, text_fontsize, fileformat


# Plot attributes
_fig_size, _title_fontsize, _subtitle_fontsize, _xy_fontsize, _xy_ticks, _legend_fontsize, _text_fontsize, _fileformat = \
    plot_properties()
