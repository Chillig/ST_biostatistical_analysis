"""Statistics: (Weighted) Pearson Correlation, Wilconxon rank sum test, Shapiro test
    File name: corr_statistics.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

import numpy as np
from scipy import stats
import scipy.special as special
from scipy.stats import rankdata
import statsmodels.api as sm

import scipy
from statsmodels.stats.multicomp import pairwise_tukeyhsd


def weighted_cov(x, y, w):
    """Weighted Covariance
        Source: https://stackoverflow.com/questions/38641691/weighted-correlation-coefficient-with-pandas

    Parameters
    ----------
    x : numpy.array, pandas.Dataframe
        responder counts
    y : numpy.array, pandas.Dataframe
        cytokine counts
    w : numpy.array, pandas.Dataframe
        size of cyto+ clusters

    Returns
    -------
    weighted covariance

    """
    return np.sum(w * (x - np.average(x, weights=w)) * (y - np.average(y, weights=w))) / np.sum(w)


def weighted_corr(x, y, w):
    """Weighted Correlation
        Source: https://stackoverflow.com/questions/38641691/weighted-correlation-coefficient-with-pandas

    Parameters
    ----------
    x :  numpy.array, pandas.Dataframe
        responder counts
    y : numpy.array, pandas.Dataframe
        cytokine counts
    w : numpy.array, pandas.Dataframe
        size of cyto+ clusters

    Returns
    -------
    weighted correlation value

    """
    return weighted_cov(x, y, w) / np.sqrt(weighted_cov(x, x, w) * weighted_cov(y, y, w))


def _wrank(x, w):
    (unique, arr_inv, counts) = np.unique(rankdata(x), return_counts=True, return_inverse=True)
    a = np.bincount(arr_inv, w)
    return (np.cumsum(a) - a)[arr_inv]+((counts + 1)/2 * (a/counts))[arr_inv]


def calculate_weighted_spearman(x, y, w):
    # taken from https://github.com/matthijsz/weightedcorr/blob/master/WeightedCorr.py
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

    # stats.probplot(df_highcounts, dist="norm", plot=plt)
    # plt.title("High Counts Q-Q Plot")
    #
    # stats.probplot(df_lowcounts, dist="norm", plot=plt)
    # plt.title("Lower Counts Q-Q Plot")

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


def apply_statstest(df_highcounts, df_lowcounts, correlation_value):
    """Check if data is normal distributed and calculate the p-value for the correlation value

    Parameters
    ----------
    df_highcounts : pandas.Dataframe
    df_lowcounts : pandas.Dataframe
    correlation_value : float

    Returns
    -------
    prob : float
        probability for the correlation value

    """
    # check if data is normal distributed
    w_hc, ps_hc, w_lc, ps_lc = check_data_normaldistributed(df_highcounts, df_lowcounts)

    # if p-value < 0.05 -> variable violates the assumption of normality => Use test
    if (ps_hc <= 0.05) & (ps_lc <= 0.05):
        ab = len(df_highcounts.values) / 2 - 1
        prob = 2 * special.btdtr(ab, ab, 0.5 * (1 - abs(np.float64(correlation_value))))
        print(prob)
        return prob
    else:
        ab = len(df_highcounts.values) / 2 - 1
        prob = 2 * special.btdtr(ab, ab, 0.5 * (1 - abs(np.float64(correlation_value))))
        print(prob)
        print('Causion: Not normally distributed data')
        return prob


def create_ls_model(x, y, w=None, const=False):
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


def calculate_errorinterval(func, x_values, y_values=None, cov=None, popt=None):
    # Confidence interval
    sigma_ab = np.sqrt(np.diagonal(cov))

    # prepare confidence level curves
    nstd = 5.  # to draw 5-sigma intervals
    popt_up = popt + nstd * sigma_ab
    popt_dw = popt - nstd * sigma_ab

    if type(func) == 'function':
        fit = func(x_values, *popt)
        fit_up = func(x_values, *popt_up)
        fit_dw = func(x_values, *popt_dw)
    else:
        stdev = np.sqrt(sum((func.predict(x_values) - y_values) ** 2) / (len(y_values) - 2))
        prediction = func.predict(x_values)
        fit_dw = prediction - 1.96 * stdev
        fit_up = prediction + 1.96 * stdev

    return fit, fit_up, fit_dw


def get_1sigma_confidenceinterval(func, x_values, cov, popt):
    # Confidence interval
    sigma_ab = np.sqrt(np.diagonal(cov))
    bound_upper = func(x_values, *(popt + sigma_ab))
    bound_lower = func(x_values, *(popt - sigma_ab))

    return bound_upper, bound_lower


def get_correlation_stats(df, resp_name, cyto, weight_name, dict_corr):
    for method in ['pearson', 'spearman']:
        if method == 'pearson':
            # Weighted Pearson Correlation (WPC)
            correlation = weighted_corr(
                x=df[resp_name], y=df[cyto], w=df[weight_name])
        else:
            # Weighted Spearman Correlation (WSC)
            correlation = calculate_weighted_spearman(
                x=df[resp_name], y=df[cyto], w=df[weight_name])
        # 1.3 Calculate p-value
        p_w = apply_statstest(
            df_highcounts=df[resp_name], df_lowcounts=df[cyto], correlation_value=correlation)

        # Save correlation value and p-value in dictionary
        dict_corr[method].append((cyto, resp_name, correlation, p_w))

    return dict_corr


def apply_hfdr(df, df_melted):
    groups_pvalues = {}
    # Check if data is normal distributed
    _, pval_ifng = scipy.stats.shapiro(df['IFNG_cluster'][~df['IFNG_cluster'].isna()])
    _, pval_il13 = scipy.stats.shapiro(df['IL13_cluster'][~df['IL13_cluster'].isna()])
    _, pval_il17a = scipy.stats.shapiro(
        df['IL17A_cluster'][~df['IL17A_cluster'].isna()])
    if (pval_ifng < 0.05) and (pval_il13 < 0.05) and (pval_il17a < 0.05):
        # Apply anova
        _, pval = scipy.stats.f_oneway(
            df['IFNG_cluster'][~df['IFNG_cluster'].isna()],
            df['IL13_cluster'][~df['IL13_cluster'].isna()],
            df['IL17A_cluster'][~df['IL17A_cluster'].isna()])
        if pval < 0.05:
            # Perform Tukey test
            tukey = pairwise_tukeyhsd(endog=df_melted['value'].astype(float), groups=df_melted['variable'], alpha=0.05)
    else:
        # Test if groups differ significantly scipy.stats.kruskal
        _, pval = scipy.stats.kruskal(
            df['IFNG_cluster'][~df['IFNG_cluster'].isna()],
            df['IL13_cluster'][~df['IL13_cluster'].isna()],
            df['IL17A_cluster'][~df['IL17A_cluster'].isna()])
        # If p-value < 0.05 -> Null hypothesis rejected => groups differ

    # TODO add tests and store p-values

    return groups_pvalues
