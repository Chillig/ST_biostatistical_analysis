import scanpy as sc


def scaling(adata):
    """Scale data to unit variance and zero mean.
        “uninteresting” sources of variation: technical noise, batch effects,
        and/or uncontrolled biological variation (e.g. cell cycle)

    Parameters
    ----------
    adata : annData

    Returns
    -------
    adata : annData

    """
    adata = sc.pp.scale(adata, zero_center=True, copy=True)

    return adata


def apply_regression_variables(adata):
    """Apply regression
    Regressing variation due to uninteresting sources can improve downstream identification
    of principal components and clustering.

    ATTENTION: Apply only if you want to neglect cell cycle biological covariates and
               don't want to look at the global structure of your data set
               see: https://github.com/theislab/scanpy/issues/722

    Parameters
    ----------
    adata : annData

    Returns
    -------
    adata : annData

    """
    # Keys for observation annotation on which to regress on.
    adata = sc.pp.regress_out(adata=adata, keys='phase', copy=True)

    return adata
