import os
import pandas as pd
import numpy as np

from scanpy import logging as logg
from scanpy.neighbors import neighbors
from anndata import AnnData
from scipy.sparse import csc_matrix

import scrublet as scr
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib


def scrublet_algorithm(adata, sample_name, save_folder):
    """Apply Scrublet algorithm to detect doublets
    1. When working with data from multiple samples, run Scrublet on each sample separately.
    Because Scrublet is designed to detect technical doublets formed by the random co-encapsulation of two cells,
    it may perform poorly on merged datasets where the cell type proportions are not representative of any single sample.

    2. Check that the doublet score threshold is reasonable
    (in an ideal case, separating the two peaks of a bimodal simulated doublet score histogram, as in this example),
    and adjust manually if necessary.

    3.Visualize the doublet predictions in a 2-D embedding (e.g., UMAP or t-SNE).
    Predicted doublets should mostly co-localize
    (possibly in multiple clusters). If they do not, you may need to adjust the doublet score threshold,
    or change the pre-processing parameters to better resolve the cell states present in your data.

    Parameters
    ----------
    adata : annData
    sample_name : str
    save_folder : str

    Returns
    -------
    adata : annData

    """
    matplotlib.use('Agg')

    cp_adata = adata.copy()
    # todo optimize simulation
    # sim_doublet_ratio: 2, n_neighbors of KNN: default 0.5 * sqrt(n_cells)

    scrub = scr.Scrublet(cp_adata.X, expected_doublet_rate=0.1)
    # scrub = scr.Scrublet(cp_adata.X)

    # Run the default pipeline, which includes:
    # 1. Doublet simulation
    # 2. Normalization, gene filtering, rescaling, PCA
    # 3. Doublet score calculation
    # 4. Doublet score threshold detection and doublet calling

    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, use_approx_neighbors=True,
                                                              min_gene_variability_pctl=85,
                                                              n_prin_comps=30, log_transform=True,
                                                              get_doublet_neighbor_parents=True, mean_center=True)

    # The parent transcriptomes that generated the doublet neighbors of each observed transcriptome. Information can
    # be used to infer the cell states that generated a given doublet state.

    scrub.call_doublets(threshold=0.25)
    fig, ax = scrub.plot_histogram()
    fig.savefig(os.path.join(save_folder, 'obs_sim_doublets.png'), bbox_inches='tight')

    # Comparison of scrublet and scanpy plotting and embeddings (the plots should look the same)
    sc.pp.highly_variable_genes(cp_adata, flavor='cell_ranger', n_top_genes=4000, batch_key='sample')
    sc.pp.pca(cp_adata, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(cp_adata, n_neighbors=10, metric='euclidean')
    sc.tl.umap(cp_adata, min_dist=0.5,  n_components=2)
    sc.pl.umap(cp_adata, color='sample')
    plt.savefig(os.path.join(save_folder, 'scanpy_UMAP.png'), bbox_inches='tight')

    # plot count matrix
    scrub.set_embedding('UMAP', scr.get_umap(cp_adata.X, n_neighbors=10, min_dist=0.5))
    fig_umap, ax_umap = scrub.plot_embedding('UMAP', order_points=True, score='zscore')
    fig_umap.savefig(os.path.join(save_folder, "_".join([sample_name, 'scrublet_zscore_count_matrix.png'])),
                     bbox_inches='tight')

    scrub.set_embedding('UMAP', scr.get_umap(cp_adata.X, n_neighbors=10, min_dist=0.5))
    fig_umap, ax_umap = scrub.plot_embedding('UMAP', order_points=True)  # , score='zscore')
    fig_umap.savefig(os.path.join(save_folder, "_".join([sample_name, 'scrublet_count_matrix.png'])),
                     bbox_inches='tight')

    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.5))
    fig_umap, ax_umap = scrub.plot_embedding('UMAP', order_points=True,)  # score='zscore')
    fig_umap.savefig(os.path.join(save_folder, "_".join([sample_name, 'scrublet_UMAP_1.png'])), bbox_inches='tight')

    doublets_barcode = np.asarray(cp_adata.obs[predicted_doublets].index)
    if len(doublets_barcode) > 0:
        doublet_df = pd.DataFrame({'doublet_barcodes': doublets_barcode, 'score': doublet_scores[predicted_doublets]})
        doublet_df.to_csv(os.path.join(save_folder, 'HVG_doublet_detection.csv'))

    # return annData object cleaned from doublets :)
    adata = adata[~predicted_doublets]

    matplotlib.use('MacOSX')  # not an ideal solution as it is system dependend ..

    return adata


# Carlos doublet detection
def woublet(adata, sim_doublet_ratio=2, n_neighbors=30, expected_doublet_rate=0.1, total_counts_key='n_counts', copy=False):
    adata = adata.copy() if copy else adata

    if 'X_pca' not in adata.obsm_keys():
        raise ValueError(
            'Did not find \'X_pca\'. Run `sc.pp.pca` first.')

    if total_counts_key in adata.obs:
        total_counts = np.array(adata.obs[total_counts_key])
    else:
        total_counts = np.ones(adata.X.shape[0])

    # Simulate doublets by averaging PC coordinates of random cell pairs
    logg.info('Simulating doublets')
    PCdat, doub_labels, parent_ix = simulate_doublets_from_pca(adata.obsm['X_pca'],
                                                               total_counts=total_counts,
                                                               sim_doublet_ratio=sim_doublet_ratio)

    adata_doub = AnnData(csc_matrix((PCdat.shape[0], 1)))
    adata_doub.obsm['X_pca'] = PCdat

    # Calculate doublet scores using k-nearest-neighbor classifier
    logg.info('Running KNN classifier')
    adata.obs['doublet_score'], adata.uns['sim_doublet_score'] = calculate_doublet_scores(
        adata_doub,
        doub_labels,
        n_neighbors=n_neighbors,
        expected_doublet_rate=expected_doublet_rate)

    # adata.obs['doublet_score']

    return adata if copy else None


# ========================================================================================#

def simulate_doublets_from_pca(PCdat, total_counts=None, sim_doublet_ratio=1):
    """
    Simulate doublets by averaging PCA coordinates of random cell pairs.
    Average is weighted by total counts of each parent cell, if provided.

    Returns:
    PCdoub (matrix of size (num_cells+num_sim_doubs, num_pcs)): PCA matrix with the simulated doublet
    PCA coordinates appended to the original data matrix PCdat.
    doub_labels (array of size (num_cells+num_sim_doubs)): 0 if observed cell, 1 if simulated doublet
    pair_ix (matrix of size(num_sim_doubs, 2)): each row gives the indices of the parent cells used to generate
    the simulated doublet
    """

    n_obs = PCdat.shape[0]
    n_doub = int(n_obs * sim_doublet_ratio)

    if len(total_counts) == 0:
        total_counts = np.ones(n_obs)

    pair_ix = np.random.randint(0, n_obs, size=(n_doub, 2))

    pair_tots = np.hstack((total_counts[pair_ix[:, 0]][:, None], total_counts[pair_ix[:, 1]][:, None]))
    pair_tots = np.array(pair_tots, dtype=float)
    pair_fracs = pair_tots / np.sum(pair_tots, axis=1)[:, None]

    PCdoub = PCdat[pair_ix[:, 0], :] * pair_fracs[:, 0][:, None] + PCdat[pair_ix[:, 1], :] * pair_fracs[:, 1][:, None]

    PCdoub = np.vstack((PCdat, PCdoub))
    doub_labels = np.concatenate((np.zeros(n_obs), np.ones(n_doub)))

    return PCdoub, doub_labels, pair_ix


# ========================================================================================#

def calculate_doublet_scores(adata, doub_labels, n_neighbors=30, expected_doublet_rate=1.0):
    n_obs = sum(doub_labels == 0)
    n_sim = sum(doub_labels == 1)

    # Adjust k (number of nearest neighbors) based on the ratio of simulated to observed cells
    k_adj = int(round(n_neighbors * (1 + n_sim / float(n_obs))))

    # Find k_adj nearest neighbors
    neighbors(adata, n_neighbors=k_adj, use_rep='X_pca')

    # Calculate doublet score based on ratio of simulated cell neighbors vs. observed cell neighbors
    matrix = adata.uns['neighbors']['distances']
    n_sim_neigh = (matrix[:, doub_labels == 1] > 0).sum(1).A.squeeze()
    n_obs_neigh = (matrix[:, doub_labels == 0] > 0).sum(1).A.squeeze()

    doub_score = n_sim_neigh / (n_sim_neigh + n_obs_neigh * n_sim / float(n_obs) / expected_doublet_rate)

    # return doublet scores for observed cells and simulated cells
    return doub_score[doub_labels == 0], doub_score[doub_labels == 1]


def carlos_woublet(adata, save_folder, key):
    """Perform scrublet like Carlos -> Woublet :)

    Parameters
    ----------
    adata : annData
    save_folder : str
    key : str

    Returns
    -------

    """
    cp_adata = adata.copy()
    holder = np.zeros((cp_adata.shape[0],))
    for smp in np.unique(cp_adata.obs['sample']):
        adata_smp = cp_adata[cp_adata.obs['sample'] == smp]
        sc.tl.pca(adata_smp, svd_solver='arpack')
        woublet(adata_smp)
        holder[cp_adata.obs['sample'] == smp] = adata_smp.obs['doublet_score']

    cp_adata.obs['doublet_score'] = holder

    sc.tl.pca(cp_adata, svd_solver='arpack')
    sc.pp.neighbors(cp_adata, method='umap', metric='euclidean')
    sc.tl.umap(cp_adata, init_pos='spectral')

    sc.settings.set_figure_params(dpi=160, dpi_save=300, vector_friendly=True, format='svg')
    cp_adata.obs['doublet_score_0.65'] = cp_adata.obs['doublet_score'] > 0.65
    cp_adata.obs['doublet_score_0.35'] = cp_adata.obs['doublet_score'] > 0.35
    cp_adata.obs['doublet_score_0.20'] = cp_adata.obs['doublet_score'] > 0.20

    sc.pl.umap(cp_adata, color=['doublet_score', 'doublet_score_0.65', 'doublet_score_0.35'],
               size=1, legend_fontsize=6)
    plt.savefig(os.path.join(save_folder, '{}_Woublet_DD.png'.format(key)), bbox_inches='tight')
    plt.close()

    print("% of cells with a doublet score above 0.65: ", cp_adata.obs['doublet_score_0.65'].sum() / cp_adata.shape[0])
    print("% of cells with a doublet score above 0.35: ", cp_adata.obs['doublet_score_0.35'].sum() / cp_adata.shape[0])
    print("% of cells with a doublet score above 0.20: ", cp_adata.obs['doublet_score_0.20'].sum() / cp_adata.shape[0])

    sc.pl.umap(cp_adata, color=['doublet_score', 'doublet_score_0.35'], size=1, frameon=False)
    plt.savefig(os.path.join(save_folder, '{}_Woublet_DD_0_35.png'.format(key)), bbox_inches='tight')
    plt.close()

    data_raw = adata[cp_adata.obs_names]
    data_raw.obs['scrublet_score'] = cp_adata.obs['doublet_score']
    print(data_raw)

    return data_raw
