import scanpy as sc
import os
import numpy as np
from datetime import date
import matplotlib.pyplot as plt
from matplotlib import colors

import scipy.stats as stats


from python_scripts.utils import gene_lists
from python_scripts.spatial_correlation import tools


def calculate_ncounts_ngenes(adata):
    """

    Parameters
    ----------
    adata : annData

    Returns
    -------
    adata : annData
        contains now observables n_counts, log_counts, n_genes and mt_frac

    """
    # sum up the number of UMI counts (barcodes) in count matrix .X along cols (genes)
    # number of counts per spot
    adata.obs['n_counts'] = adata.X.sum(1)
    # number of counts per spot in log
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    # number of genes per spot
    adata.obs['n_genes'] = (adata.X > 0).sum(1)

    mt_gene_mask = [gene.startswith('MT-') for gene in adata.var_names]
    adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1) / adata.obs['n_counts']

    # add ribosomal genes
    # mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    # ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var['hb'] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], percent_top=None, log1p=True, inplace=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], percent_top=None, log1p=False, inplace=True)

    print("Overview of n_counts: ", adata.obs['n_counts'].describe())
    print("Overview of n_genes: ", adata.obs['n_genes'].describe())

    return adata


def plot_scatter(adata, title, cyto_resps, resp_label, save_folder):
    fig, ax = plt.subplots(figsize=(12, 6), ncols=2, nrows=1)
    sa = ax[0].scatter(adata.obs["{}_counts".format(resp_label)], adata.obs['mt_frac'], c=adata.obs['n_counts'], s=3)
    ax[0].set_xlabel("Responder genes counts")
    ax[0].set_ylabel("MT-fraction")
    ax[0].spines["top"].set_visible(False)
    ax[0].spines["right"].set_visible(False)
    cb = fig.colorbar(sa, ax=ax[0])
    cb.set_label('Number of counts')
    ax[0].set_title(cyto_resps)

    sa = ax[1].scatter(adata.obs["{}_counts".format(resp_label)], adata.obs['mt_frac'], c=adata.obs['n_genes'], s=3)
    ax[1].set_xlabel("Responder genes counts")
    ax[1].set_ylabel("MT-fraction")
    ax[1].spines["top"].set_visible(False)
    ax[1].spines["right"].set_visible(False)
    cb = fig.colorbar(sa, ax=ax[1])
    cb.set_label('Number of genes')
    ax[1].set_title(cyto_resps)

    plt.tight_layout()

    plt.savefig(os.path.join(save_folder, "{}__{}_MT_vs_responder_counts.pdf".format(cyto_resps, title)))
    plt.close()


def plot_mt_counts_gene_filter(adata_responders, cytokine, responders, save_folder):
    fig, ax = plt.subplots(figsize=(12, 6), ncols=2, nrows=1)
    # min UMI-counts = 50
    sa = ax[0].scatter(adata_responders.obs["{}_counts".format(responders)], adata_responders.obs['mt_frac'],
                       c=adata_responders.obs['n_counts'], s=4)
    mask = (adata_responders.obs['n_counts'] > 50) & (adata_responders.obs['mt_frac'] > 0.25)
    ax[0].scatter(adata_responders[mask].obs["{}_counts".format(responders)],
                  adata_responders[mask].obs['mt_frac'], c='red', s=3)
    mask = adata_responders.obs['n_counts'] < 50
    ax[0].scatter(adata_responders[mask].obs["{}_counts".format(responders)],
                  adata_responders[mask].obs['mt_frac'], c='grey', s=6)

    ax[0].axhline(y=0.25, color='darkred', linestyle='--')
    ax[0].set_xlabel("Responder genes counts")
    ax[0].set_ylabel("MT-fraction")
    ax[0].spines["top"].set_visible(False)
    ax[0].spines["right"].set_visible(False)

    cb = fig.colorbar(sa, ax=ax[0])
    cb.set_label('Number of counts')
    ax[0].set_title(cytokine)

    # min genes = 30
    sa = ax[1].scatter(adata_responders.obs["{}_counts".format(responders)], adata_responders.obs['mt_frac'],
                       c=adata_responders.obs['n_genes'], s=4)
    mask = (adata_responders.obs['n_genes'] > 30) & (adata_responders.obs['mt_frac'] > 0.25)
    ax[1].scatter(adata_responders[mask].obs["{}_counts".format(responders)],
                  adata_responders[mask].obs['mt_frac'], c='red', s=3)
    mask = adata_responders.obs['n_genes'] < 30
    ax[1].scatter(adata_responders[mask].obs["{}_counts".format(responders)],
                  adata_responders[mask].obs['mt_frac'], c='grey', s=6)
    ax[1].axhline(y=0.25, color='darkred', linestyle='--')
    ax[1].set_xlabel("Responder genes counts")
    ax[1].set_ylabel("MT-fraction")
    ax[1].spines["top"].set_visible(False)
    ax[1].spines["right"].set_visible(False)
    cb = fig.colorbar(sa, ax=ax[1])
    cb.set_label('Number of genes')
    ax[1].set_title(cytokine)

    plt.tight_layout()

    plt.savefig(os.path.join(save_folder, "{}__MT_vs_responder_counts_n_countsngenes.pdf".format(cytokine)))
    plt.close()


def density_estimation(m1, m2):
    xmin = np.amin(m1)
    xmax = np.amax(m1)
    ymin = np.amin(m2)
    ymax = np.amax(m2)
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z


def main(save_folder):
    adata = sc.read('/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/adata_storage/2021-07-29/Spatial Transcriptomics_unpp.h5')

    adata = calculate_ncounts_ngenes(adata=adata)

    # remove zero count spots ..
    adata = adata[adata.obs['n_counts'] != 0]

    # Inflammatory markers == responder genes of cytokines
    _, _, cyto_responders = gene_lists.get_publication_cyto_resps()

    # counts of cytokine (cyto_resps)
    adata = tools.add_columns_genes(adata=adata, genes=list(cyto_responders.keys()), label='cytokines')

    mask_allcuts = (adata.obs['n_counts'] < 50) | (adata.obs['mt_frac'] > 0.25) | (adata.obs['n_genes'] < 30)
    adata_cutted = adata[mask_allcuts].copy()
    number_cytokine_removed = len(np.where(adata_cutted.obs['cytokines_label'] == 'cytokines')[0])
    print("Number of Cytokine positive spots which get filtered out by using our Threshold "
          "n_counts > 50, MT-threshold < 25%, and n_genes > 30: ", number_cytokine_removed)

    fig, ax = plt.subplots(figsize=(12, 6), ncols=2, nrows=1)
    # min UMI-counts = 50
    ax[0].scatter(adata.obs['n_counts'], adata.obs['mt_frac'], s=6, alpha=0.3, edgecolors='none', label='pass cut-offs')
    # # add contour lines
    # X, Y, Z = density_estimation(adata.obs['n_counts'], adata.obs['mt_frac'])
    # ax[0].contour(X, Y, Z, 6, colors='darkgrey')

    # spots which are removed due to the MT-threshold
    mask = (adata.obs['n_counts'] > 50) & (adata.obs['mt_frac'] > 0.25)
    ax[0].scatter(adata[mask].obs['n_counts'], adata[mask].obs['mt_frac'], c='red', s=6, alpha=0.3, edgecolors='none',
                  label='removed by MT-fraction cut-off > 25%')
    ax[0].text(250000, 0.26, 'Number of spots: {}'.format(adata[mask].shape[0]))
    ax[0].text(250000, 0.22, 'Number of spots: {}'.format(adata[~mask].shape[0]))
    # spots which are filtered out by n_counts
    mask = adata.obs['n_counts'] < 50
    ax[0].scatter(adata[mask].obs['n_counts'], adata[mask].obs['mt_frac'], c='grey', s=6, alpha=0.3, edgecolors='none',
                  label='filtered by n_counts < 50')
    # Highlight cyto+ spots
    mask = adata.obs['cytokines_counts'] > 0
    ax[0].scatter(adata[mask].obs['n_counts'], adata[mask].obs['mt_frac'], c='darkorange', s=6, alpha=1,
                  edgecolors='none', label='cytokines')
    # Highlight cyto+ spots which would be removed by filtering n_counts
    mask_cyto_removed = (adata.obs['cytokines_counts'] > 0) & (adata.obs['n_counts'] < 50)
    ax[0].scatter(adata[mask_cyto_removed].obs['n_counts'], adata[mask_cyto_removed].obs['mt_frac'],
                  s=4, alpha=1, edgecolors='grey', facecolors='none', linewidth=1)

    ax[0].axhline(y=0.25, color='darkred', linestyle='--')
    ax[0].axvline(x=50, color='darkred', linestyle='--')
    ax[0].set_xlabel("Number of UMI-counts")
    ax[0].set_ylabel("MT-fraction")
    ax[0].spines["top"].set_visible(False)
    ax[0].spines["right"].set_visible(False)
    ax[0].set_title('Min UMI-counts filter 50 and max. MT-fraction cut-off of 25 %')
    ax[0].legend(frameon=False)

    # min genes = 30
    ax[1].scatter(adata.obs['n_genes'], adata.obs['mt_frac'], s=6, alpha=0.3, edgecolors='none', label='pass cut-offs')
    # spots which are removed due to the MT-threshold
    mask = (adata.obs['n_genes'] > 30) & (adata.obs['mt_frac'] > 0.25)
    ax[1].scatter(adata[mask].obs['n_genes'], adata[mask].obs['mt_frac'], c='red', s=6, alpha=0.3, edgecolors='none',
                  label='removed by MT-fraction cut-off > 25%')
    ax[1].text(7000, 0.26, 'Number of spots: {}'.format(adata[mask].shape[0]))
    ax[1].text(7000, 0.22, 'Number of spots: {}'.format(adata[~mask].shape[0]))
    # spots which are filtered out by n_genes
    mask = adata.obs['n_genes'] < 30
    ax[1].scatter(adata[mask].obs['n_genes'], adata[mask].obs['mt_frac'], c='grey', s=4, alpha=0.3, edgecolors='none',
                  label='filtered by n_genes < 30')
    # Highlight cyto+ spots
    mask_cyto = adata.obs['cytokines_counts'] > 0
    ax[1].scatter(adata[mask_cyto].obs['n_genes'], adata[mask_cyto].obs['mt_frac'], c='darkorange', s=6, alpha=1,
                  edgecolors='none', label='cytokines')
    # Highlight cyto+ spots which would be removed
    mask_cyto_removed = (adata.obs['cytokines_counts'] > 0) & (adata.obs['n_genes'] < 30)
    ax[1].scatter(adata[mask_cyto_removed].obs['n_genes'], adata[mask_cyto_removed].obs['mt_frac'],
                  s=4, alpha=1, edgecolors='grey', facecolors='none', linewidth=1)
    ax[1].axhline(y=0.25, color='darkred', linestyle='--')
    ax[1].axvline(x=30, color='darkred', linestyle='--')
    ax[1].set_xlabel("Number of genes")
    ax[1].set_ylabel("MT-fraction")
    ax[1].spines["top"].set_visible(False)
    ax[1].spines["right"].set_visible(False)
    ax[1].set_title('Min genes filter 30 and max. MT-fraction cut-off of 25 %')
    ax[1].legend(frameon=False)

    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "MT-fraction_desicion.pdf"))
    plt.close()


    # colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
    # colors3 = plt.cm.Greys_r(np.linspace(0.7, 0.8, 1))
    # colorsComb = np.vstack([colors3, colors2])
    # mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

    # get counts of responder genes per cytokine
    resp_label = []
    for ind, cyto_resps in enumerate(cyto_responders.keys()):
        # counts of signature responders of a cytokine
        resp_label.append("_".join([cyto_resps, "Responders"]))
        adata = tools.add_columns_genes(adata=adata, genes=cyto_responders[cyto_resps],
                                        label=resp_label[ind])

        adata_responders = adata[adata.obs["{}_counts".format(resp_label[ind])] > 0].copy()

        # plot scatter plot MT-fraction vs. responder genes counts
        plot_mt_counts_gene_filter(
            adata_responders=adata_responders, responders=resp_label[ind], save_folder=save_folder, cytokine=cyto_resps)

        # Plot those spots which would be filtered out by QC-Filtering due to MT-fraction
        mask = \
            ((adata_responders.obs['n_genes'] > 30) & (adata_responders.obs['n_counts'] > 50)) & \
            (adata_responders.obs['mt_frac'] > 0.25)

        adata_qc = adata_responders[mask].copy()
        plot_scatter(adata=adata_qc, title='Filtered_by_MT', cyto_resps=cyto_resps, resp_label=resp_label[ind],
                     save_folder=save_folder)

        # Plot those spots which would be filtered out by QC-Filtering due to n_counts or n_genes
        # but which also have a MT-fraction > 25 %
        mask = \
            ((adata_responders.obs['n_genes'] < 30) | (adata_responders.obs['n_counts'] < 50)) & \
            (adata_responders.obs['mt_frac'] > 0.25)

        adata_qc = adata_responders[mask].copy()
        plot_scatter(adata=adata_qc, title='Filtered_by_ngenes_ncounts', cyto_resps=cyto_resps,
                     resp_label=resp_label[ind], save_folder=save_folder)


if __name__ == '__main__':
    # TOdo show MT-fraction distriubution per tissue slide on H&E image
    today = date.today()
    path = os.path.join("..", "..")
    # create saving folder in current project path
    savepath = os.path.join(path, "output", "MT-fraction", str(today))
    os.makedirs(savepath, exist_ok=True)

    main(save_folder=savepath)
