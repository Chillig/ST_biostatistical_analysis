from python_scripts.pre_processing.plot_functions import plots_preprocessing

import scanpy as sc


def calc_visualization(adata, save_folder, batch_key, n_comps=50, use_highly_variable=True,
                       svd_solver='arpack', n_jobs=12, raw=False):
    """Dimension reduce, compute neighorhood graph amd embedded data

    Parameters
    ----------
    adata : annData
    save_folder : str
    batch_key : str
    n_comps : int
    use_highly_variable : bool
    svd_solver : str
    n_jobs : int
    raw : bool

    Returns
    -------

    """
    # dimension reduction
    sc.pp.pca(adata, n_comps=n_comps, use_highly_variable=use_highly_variable, svd_solver=svd_solver)

    # Compute a neighborhood graph of observations
    # metric can be: * euclidean, * manhattan, * chebyshev, * minkowski, * canberra, * braycurtis,
    # * mahalanobis, * wminkowski, * seuclidean, * cosine, * correlation, * haversine, * hamming, * jaccard, * dice,
    # * russelrao, * kulsinski, * rogerstanimoto, * sokalmichener, * sokalsneath, * yule
    sc.pp.neighbors(adata, method='umap', metric='euclidean')

    # embeddings
    # Note n_jobs works for MulticoreTSNE, but not regular implementation)
    sc.tl.tsne(adata, n_jobs=n_jobs, perplexity=30)  # t-sne very sensitive to perplexity value [5, 50]
    sc.tl.umap(adata, init_pos='spectral')
    sc.tl.diffmap(adata, neighbors_key='neighbors')
    sc.tl.draw_graph(adata, layout='fa')

    plots_preprocessing.plot_visualization_results(adata=adata, save_folder=save_folder, batch_key=batch_key, raw=raw)

    return adata
