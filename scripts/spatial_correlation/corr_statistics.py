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


def check_data_normaldist(df_highcounts, df_lowcounts):
    """Check with Shapiro test if count data is normal distributed

    Parameters
    ----------
    df_highcounts : pandas.Dataframe
    df_lowcounts : pandas.Dataframe

    Returns
    -------
    ps_hc : float
        p-value for high counts distribution
    ps_lc: float
        p-value for low counts distribution

    """
    # check if data is normal distributed
    w, ps_hc = stats.shapiro(df_highcounts)
    w, ps_lc = stats.shapiro(df_lowcounts)
    return ps_hc, ps_lc


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
    ps_hc, ps_lc = check_data_normaldist(df_highcounts, df_lowcounts)

    # if p-value < 0.05 -> variable violates the assumption of normality => Use Wilcoxon signed rank test
    if (ps_hc <= 0.05) & (ps_lc <= 0.05):
        t, prob = stats.wilcoxon(df_highcounts, df_lowcounts)
        print(prob)
        return prob
    else:
        return None


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
    ps_hc, ps_lc = check_data_normaldist(df_highcounts, df_lowcounts)

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
