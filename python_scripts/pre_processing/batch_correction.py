import python_scripts.pre_processing.plot_functions.plots_preprocessing as plt_pp_plots
import python_scripts.pre_processing.pc_determination as pc_determination

# Package for batch correction
import scanorama

# Operation and analysis packages
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ann

# ML models
from sklearn.linear_model import LinearRegression
import scipy.sparse as sparse

"""
Apply batch correction to correct for technical variation which can occur in the lab
Benchmark study on batch correction methods
https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-019-1850-9
results:
Best batch correction methods: LIGER, Harmony, and Seurat 3.
Harmony performed well on datasets with common cell types, and also different technologies.
LIGER performed well, especially on datasets with non-identical cell types.
To improve recovery of DEGs in batch-corrected data, use scMerge for batch correction. (Lotfollahi et al., 2019)s
Evaluate batch correction using the four assessment metrics: ASW, ARI, LISI, and kBET
"""


def _score_uncorrectd_data(adatas, n_comps, save_folder, possible_batch_effects, batch_key='library_id'):
    """

    Parameters
    ----------
    adatas : annData
    n_comps : int
    save_folder : str
    possible_batch_effects : list, str
    batch_key : str

    Returns
    -------

    """
    # Score uncorrectd matrix
    pca = pc_determination.pcs_combs(adatas.X, save_folder, raw=False, type_dataset="No_HVG_uncorrected",
                                     use_highly_variable=False, copy=True, return_info=True)

    # Score embedding
    dict_r2var = dict()
    # Score variance contribution by batch_key
    score_variance(adata=adatas, obs=batch_key, pca=pca, r2var=dict_r2var)
    for poss_be in possible_batch_effects:
        # Score variance contribution by other covariate
        score_variance(adata=adatas, obs=poss_be, pca=pca, r2var=dict_r2var)
    print(dict_r2var)

    # Plot uncorrected matrix
    sc.pp.pca(adatas, n_comps=n_comps, use_highly_variable=False, svd_solver='arpack')
    sc.pp.neighbors(adatas)
    sc.tl.umap(adatas)

    # Plot
    for covar in possible_batch_effects:
        plt_pp_plots.plot_batch_correction(adatas, save_folder, batch_key="unbc_matrix", possible_batch_effect=covar)


def _split_batches(adata, batch, hvg=None):
    """Split adata by batches
    From scIB

    Parameters
    ----------
    adata : annData
    batch : str
    hvg : None, pandas.Series

    Returns
    -------

    """

    split = []
    if hvg is not None:
        adata = adata[:, hvg]
    for i in adata.obs[batch].unique():
        split.append(adata[adata.obs[batch] == i].copy())
    return split


def score_variance(adata, obs, pca, r2var):
    """Calculate the variance of the co-variate

    Parameters
    ----------
    adata : annData
    obs : str
    pca : numpy.array
    r2var : pandas.Dataframe

    Returns
    -------

    """
    x_pca = pca[0].copy()
    n_comps = x_pca.shape[1]
    pca_sd = pca[3].copy()

    batch = adata.obs[obs].astype('category').cat.codes
    batch = pd.get_dummies(batch) if 'category' == str(batch.dtype) else np.array(batch)

    r2 = []
    for i in range(n_comps):
        lm = LinearRegression()
        lm.fit(x_pca[:, [i]], batch)
        r2.append(lm.score(x_pca[:, [i]], batch))

    var = pca_sd ** 2 / sum(pca_sd ** 2) * 100
    r2var[obs] = sum(r2 * var) / 100

    return r2var


def scanorama_bc(adatas, n_comps, save_folder, possible_batch_effects, batch_key='library_id'):
    """Apply Scanorama Batch correction
    Scanorama enables batch-correction and integration of heterogeneous scRNA-seq datasets

    Parameters
    ----------
    adatas : annData
    n_comps : int
    save_folder : str
    possible_batch_effects : list, str
    batch_key : str

    Returns
    -------

    """
    # 1. Score uncorrectd matrix
    _score_uncorrectd_data(adatas=adatas, n_comps=n_comps, save_folder=save_folder, batch_key=batch_key,
                           possible_batch_effects=possible_batch_effects)

    # 2. Apply scanoroma batch correction method
    # 2.1 Split list into annData objects
    split = _split_batches(adatas[:, adatas.var['highly_variable']].copy(), batch=batch_key)

    # 2.2 run scanorama batch correction
    kwargs = {"return_dimred": True}
    emb, corrected = scanorama.correct_scanpy(split, **kwargs)
    # concatenate corrected adatas
    emb = np.concatenate(emb, axis=0)
    adata_cor = ann.AnnData.concatenate(*corrected, batch_key=batch_key,
                                        batch_categories=adatas.obs[batch_key].cat.categories,).copy()
    adatas.obsm['X_emb'] = emb

    # 2.3 Score correct matrix
    # 2.3.2 Determine No. PCs
    pca = pc_determination.pcs_combs(adatas.obsm['X_emb'], save_folder, raw=False, type_dataset="No_HVG_corrected",
                                     use_highly_variable=False, copy=True, return_info=True)

    # 2.4 Calculate variance after batch correction - might be that the variance increased within a covariate
    dict_r2var = dict()
    adata_cor.obs[batch_key] = adatas.obs[batch_key].values
    # Score variance contribution by batch
    score_variance(adata=adata_cor, obs=batch_key, pca=pca, r2var=dict_r2var)

    for poss_be in possible_batch_effects:
        adata_cor.obs[poss_be] = adatas.obs[poss_be].values

        # Score variance contribution by other covariate
        score_variance(adata=adata_cor, obs=poss_be, pca=pca, r2var=dict_r2var)

    print(dict_r2var)

    # Compute Visualisation of corrected matrix
    try:
        n_comps = int(input("Please provide the No. principal components (default 50): "))
    except ValueError:
        n_comps = int(50)
    sc.pp.pca(adata_cor, n_comps=n_comps, use_highly_variable=False, svd_solver='arpack')
    sc.pp.neighbors(adata_cor)
    sc.tl.umap(adata_cor)

    # Plot
    for poss_be in possible_batch_effects:
        plt_pp_plots.plot_batch_correction(adata_cor, save_folder, batch_key="bc_matrix", possible_batch_effect=poss_be)

    # Compute Visualisation of corrected embedding
    sc.pp.neighbors(adatas, use_rep='X_emb')
    sc.tl.umap(adatas)

    for poss_be in possible_batch_effects:
        plt_pp_plots.plot_batch_correction(adatas, save_folder, batch_key="bc_embedding", possible_batch_effect=poss_be)

    adatas.X = sparse.csr_matrix(adatas.X)

    return adatas


def apply_batch_correction(normed_scaled_adatas, save_folder, n_comps, possible_batch_effects, batch_key='library_id'):
    """Apply Batch correction using scanorama

    Causes of Batch effects:
    - Laboratory conditions
    - Choice of reagent lot or batch
    - Personnel differences
    - Time of day when the experiment was conducted
    - Atmospheric ozone levels
    - Instruments used to conduct the experiment

    Parameters
    ----------
    normed_scaled_adatas : annData
    save_folder : str
    n_comps : int
    possible_batch_effects : list, str
    batch_key : str

    Returns
    -------

    """

    # 2.3 Batch Correction (Data integration phase)
    print("         Batch Correction")
    bc_adata = scanorama_bc(adatas=normed_scaled_adatas, n_comps=n_comps, save_folder=save_folder, batch_key=batch_key,
                            possible_batch_effects=possible_batch_effects)

    return bc_adata
