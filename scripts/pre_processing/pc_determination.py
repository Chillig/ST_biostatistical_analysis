import scripts.pre_processing.plot_functions.plots_preprocessing as pp_plt

import scanpy as sc


def pcs_combs(adata, save_folder, type_dataset, use_highly_variable, copy, return_info, raw=False):
    """Determine number of principal components for batch correction

    Parameters
    ----------
    adata : annData
    save_folder : str
    type_dataset : str
    use_highly_variable : bool
    copy : bool
    return_info : bool
    raw : bool

    Returns
    -------

    """
    if return_info:
        # takes 'numpy.ndarray' object
        pca = sc.pp.pca(adata, n_comps=50, use_highly_variable=use_highly_variable, svd_solver='arpack', copy=copy,
                        return_info=return_info)
        return pca
    else:
        # needs annData object as input
        sc.tl.pca(adata, n_comps=50, use_highly_variable=use_highly_variable, svd_solver='arpack', copy=copy,
                  return_info=return_info)

        pp_plt.plot_pc_combs(adata=adata, type_dataset=type_dataset, save_folder=save_folder, raw=raw)
