"""Statistics: (Weighted) Pearson Correlation, Wilconxon rank sum test, Shapiro test
    File name: corr_statistics.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 09/01/2025
    Python Version: 3.8
"""

import numpy as np
import pandas as pd

from scipy import stats
import scipy.special as special
from scipy.stats import rankdata
import statsmodels.api as sm


def weighted_cov(x, y, w):
    """
   Compute the weighted covariance between two variables.

   Parameters
   ----------
   x : array-like
       First variable (e.g., responder counts).
   y : array-like
       Second variable (e.g., cytokine counts).
   w : array-like
       Weights corresponding to the size or importance of each observation
       (e.g., size of cyto+ clusters).

   Returns
   -------
   float
       Weighted covariance between x and y.

   Notes
   -----
   Formula:
       cov_w(x, y) = Σ w_i * (x_i - x̄_w) * (y_i - ȳ_w) / Σ w_i
   where x̄_w and ȳ_w are the weighted means of x and y.

   Source
   ------
   https://stackoverflow.com/questions/38641691/weighted-correlation-coefficient-with-pandas
   """
    return np.sum(w * (x - np.average(x, weights=w)) * (y - np.average(y, weights=w))) / np.sum(w)


def weighted_corr(x, y, w):
    """
    Compute the weighted Pearson correlation between two variables.

    Parameters
    ----------
    x : array-like
        First variable (e.g., responder counts).
    y : array-like
        Second variable (e.g., cytokine counts).
    w : array-like
        Weights corresponding to the size or importance of each observation.

    Returns
    -------
    float
        Weighted correlation coefficient between x and y.

    Notes
    -----
    Formula:
        corr_w(x, y) = cov_w(x, y) / sqrt(cov_w(x, x) * cov_w(y, y))

    Source
    ------
    https://stackoverflow.com/questions/38641691/weighted-correlation-coefficient-with-pandas
    """
    return weighted_cov(x, y, w) / np.sqrt(weighted_cov(x, x, w) * weighted_cov(y, y, w))


def _wrank(x, w):
    """
    Compute weighted ranks of a variable for use in weighted Spearman correlation.

    Parameters
    ----------
    x : array-like
        Input variable to rank.
    w : array-like
        Weights for each observation.

    Returns
    -------
    array-like
        Weighted ranks of the input variable.
    """
    (unique, arr_inv, counts) = np.unique(rankdata(x), return_counts=True, return_inverse=True)
    a = np.bincount(arr_inv, w)
    return (np.cumsum(a) - a)[arr_inv]+((counts + 1)/2 * (a/counts))[arr_inv]


def calculate_weighted_spearman(x, y, w):
    """
    Compute the weighted Spearman correlation between two variables.

    Parameters
    ----------
    x : array-like
        First variable (e.g., responder counts).
    y : array-like
        Second variable (e.g., cytokine counts).
    w : array-like
        Weights corresponding to the size or importance of each observation.

    Returns
    -------
    float
        Weighted Spearman correlation coefficient between x and y.

    Notes
    -----
    Uses the ranking method _wrank to compute weighted ranks before calculating
    the weighted correlation. This generalizes the standard Spearman correlation
    to include weights.

    Source
    ------
    https://github.com/matthijsz/weightedcorr/blob/master/WeightedCorr.py
    """
    return weighted_corr(x=_wrank(x=x, w=w), y=_wrank(x=y, w=w), w=w)


def check_data_normaldistributed(df_highcounts, df_lowcounts):
    """Check with Shapiro test if count data is normal distributed

    Parameters
    ----------
    df_highcounts : pandas.Dataframe
    df_lowcounts : pandas.Dataframe

    Returns
    -------
    w_hc :
    ps_hc : float
        p-value for high counts distribution
    w_lc :
    ps_lc: float
        p-value for low counts distribution


    """
    # check if data is normal distributed
    w_hc, ps_hc = stats.shapiro(df_highcounts)
    w_lc, ps_lc = stats.shapiro(df_lowcounts)

    return w_hc, ps_hc, w_lc, ps_lc


def apply_wilcoxontest(df_highcounts, df_lowcounts):
    """Check if data is normally distributed

    Parameters
    ----------
    df_highcounts : pandas.Dataframe
    df_lowcounts : pandas.Dataframe

    Returns
    -------

    """

    # check if data is normal distributed
    w_hc, ps_hc, w_lc, ps_lc = check_data_normaldistributed(df_highcounts, df_lowcounts)

    # if p-value < 0.05 -> variable violates the assumption of normality => Use Wilcoxon signed rank test
    if (ps_hc <= 0.05) | (ps_lc <= 0.05):
        if len(df_highcounts) == len(df_lowcounts):
            t, prob = stats.wilcoxon(df_highcounts, df_lowcounts)
        else:
            u, prob = stats.mannwhitneyu(df_highcounts, df_lowcounts)
    else:
        stat_value, p = stats.levene(df_highcounts, df_lowcounts)
        if p < 0.05:
            # Data has not equal variance -> Welch's t-test
            w, prob = stats.ttest_ind(df_highcounts, df_lowcounts, equal_var=False)
        else:
            t, prob = stats.ttest_ind(df_highcounts, df_lowcounts)

    print(prob)
    return prob


def weighted_spearman_permutation_test(x, y, w, n_permutations: int = 10000, seed: int = None):
    """
    Compute weighted Spearman correlation and its p-value via permutation test.

    Parameters
    ----------
    x : array-like
        First variable.
    y : array-like
        Second variable.
    w : array-like
        Weights for each observation.
    n_permutations : int, default=10000
        Number of permutations for null distribution.
    seed : int or None, default=None
        Random seed for reproducibility.

    Returns
    -------
    observed_corr : float
        Weighted Spearman correlation coefficient.
    p_value : float
        Two-tailed empirical p-value from permutation test.
    """
    np.random.seed(seed)

    observed_corr = calculate_weighted_spearman(x, y, w)
    permuted_corrs = []

    for _ in range(n_permutations):
        permuted_y = np.random.permutation(y)
        perm_corr = calculate_weighted_spearman(x, permuted_y, w)
        permuted_corrs.append(perm_corr)

    permuted_corrs = np.array(permuted_corrs)
    p_value = np.mean(np.abs(permuted_corrs) >= np.abs(observed_corr))

    return observed_corr, p_value


def apply_statstest(df_highcounts: pd.DataFrame, df_lowcounts: pd.DataFrame, correlation_value: float):
    """
    Compute Pearson correlation p-value given two datasets and a correlation coefficient.
    Assumes normally distributed data. If normality is violated, returns the same
    calculation but issues a warning.

    Parameters
    ----------
    df_highcounts : pandas.DataFrame
        Gene counts or measurements for the "high expression" group.
    df_lowcounts : pandas.DataFrame
        Gene counts or measurements for the "low expression" group.
    correlation_value : float
        Pearson correlation coefficient for which significance should be tested.

    Returns
    -------
    prob : float
        Two-tailed p-value for the given Pearson correlation coefficient.
    """
    w_hc, ps_hc, w_lc, ps_lc = check_data_normaldistributed(df_highcounts, df_lowcounts)

    ab = len(df_highcounts.values) / 2 - 1
    prob = 2 * special.btdtr(ab, ab, 0.5 * (1 - abs(np.float64(correlation_value))))

    if (ps_hc <= 0.05) & (ps_lc <= 0.05):
        return prob
    else:
        print('Caution: Data not normally distributed — Pearson p-value may be invalid.')
        return prob


def get_correlation_stats(df: pd.DataFrame, proximity_genes_name: str, gene, weight_name: str, dict_corr: dict):
    """
    Compute weighted Pearson and weighted Spearman correlations with p-values.

    Weighted Pearson Correlation (WPC)
    - Correlation is weighted by transcript counts.
    - p-values are computed using a method appropriate for weighted correlations.
    - Assumes normality; warnings are displayed if normality assumptions are violated.

    Weighted Spearman Correlation (WSC)
    - Non-parametric rank-based correlation, weighted by transcript counts.
    - p-values are computed via permutation testing, creating an empirical null distribution from permuted data.

    Parameters
    ----------
    df : pandas.DataFrame
        Input data containing response, gene expression, and weights.
    proximity_genes_name : str
        Column name for genes in close vicinity.
    gene : str
        Column name for gene expression.
    weight_name : str
        Column name for weights.
    dict_corr : dict
        Dictionary to store results. Must have keys 'pearson' and 'spearman'.

    Returns
    -------
    dict_corr : dict
        Updated dictionary with results appended.
    """
    # Make sure all are float values
    df[gene] = df[gene].values.astype(np.float64)
    df[proximity_genes_name] = df[proximity_genes_name].values.astype(np.float64)
    df[weight_name] = df[weight_name].values.astype(np.float64)

    # Weighted Pearson
    correlation = weighted_corr(x=df[proximity_genes_name], y=df[gene], w=df[weight_name])
    p_w = apply_statstest(df_highcounts=df[proximity_genes_name], df_lowcounts=df[gene], correlation_value=correlation)
    dict_corr['pearson'].append((gene, proximity_genes_name, correlation, p_w))

    # Weighted Spearman
    correlation, p_w = weighted_spearman_permutation_test(
        x=df[proximity_genes_name], y=df[gene], w=df[weight_name], n_permutations=10000, seed=0)
    dict_corr['spearman'].append((gene, proximity_genes_name, correlation, p_w))

    return dict_corr


def create_ls_model(x, y, w=None, const=False):
    """
    Fit a (weighted) least squares regression model and return predictions and summary statistics.

    Parameters
    ----------
    x : array-like or pandas.DataFrame
        Predictor variable(s). If multiple predictors, provide as 2D array or DataFrame.
    y : array-like or pandas.Series
        Response variable.
    w : array-like, optional (default=None)
        Weights for weighted least squares regression. If None, ordinary least squares (OLS) is used.
    const : bool, optional (default=False)
        If True, adds a constant (intercept) term to the model.

    Returns
    -------
    est_fitted : statsmodels.regression.linear_model.RegressionResultsWrapper
        Fitted regression model object with parameter estimates, residuals, and diagnostics.
    ypred : array-like
        Predicted values from the fitted model.
    df_results : pandas.DataFrame
        Summary DataFrame of predictions including:
        - mean: predicted values
        - mean_se: standard error of predicted mean
        - mean_ci_lower, mean_ci_upper: 95% confidence interval for predicted mean
        - obs_ci_lower, obs_ci_upper: 95% prediction interval for new observations

    Notes
    -----
    - If `w` is provided, a weighted least squares (WLS) regression is performed.
    - If `const=True`, a constant term is prepended to the predictors.
    - The function prints a summary table from statsmodels `fit.summary()`.
    """
    if const:
        x = sm.add_constant(x, prepend=False)
    if w is None:
        est = sm.OLS(y, x, missing='drop')
    else:
        est = sm.WLS(y, x, w)

    est_fitted = est.fit()
    print(est_fitted.summary())

    # Fit predictions and parameters
    ypred = est_fitted.predict(x)
    pred_model = est_fitted.get_prediction(x)
    df_results = pred_model.summary_frame(alpha=0.05)

    return est_fitted, ypred, df_results
