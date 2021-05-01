import os
import pandas as pd
import numpy as np

import scrublet as scr
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


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
    # sc.pp.pca(cp_adata, n_comps=30, use_highly_variable=True)
    # sc.pp.neighbors(cp_adata, n_neighbors=10, metric='euclidean')
    # sc.tl.umap(cp_adata, min_dist=0.5,  n_components=2)
    # sc.pl.umap(cp_adata, color='sample')
    # plt.savefig(os.path.join(save_folder, 'scanpy_UMAP.png'), bbox_inches='tight')

    # plot count matrix
    scrub.set_embedding('UMAP', scr.get_umap(cp_adata.X, n_neighbors=10, min_dist=0.5))
    fig_umap, ax_umap = scrub.plot_embedding('UMAP', order_points=True, )  # score='zscore')
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

    return adata
