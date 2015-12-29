#!/usr/bin/env python

import pandas as pd
from copy import deepcopy
from . import cov_utils as cutils

class ResChar(object):

def vparameters(res):
    """
    Return the estimate of the varied parameters

    Parameters
    ----------
    res : `sncosmo.utils.Result` instance, mandatory

    """
    # parameters = deepcopy(res.parameters)
    # vparam_names = res.vparam_names
    vparameters = [res.parameters[res.param_names.index(v)]
                   for v in res.vparam_names]

    return pd.Series(vparameters)

def subCovariance(covariance, paramList, array=False):
    """
    returns the covariance of a subset of parameters in a covariance dataFrame.

    Parameters
    ----------
    covariance :

    paramList :

    Returns
    -------

    """
    arr = covariance.ix[paramList, paramList]
    if array:
        return arr
    return pd.DataFrame(arr, columns=paramList, index=paramList)




def covariance(res, normalized=False):
    """
    Returns the covariance Matrix for the varied parameters as a
    `pandas.DataFrame`, so that covariance elements may be called
    by either index or parameters.

    Parameters
    ----------
    res : `sncosmo.utils.result` instance, mandatory 

    Returns
    -------
    a `pandas.DataFrame` with column names and indexes given by the parameter
    names


    Examples
    --------
    >>> cov = covariance(res)
    >>> cov.ix[['t0', 'x1'],['t0', 'x1']]
    >>> cov.iloc[[0, 2], [0, 2]]
    """
    _cov = deepcopy(res[0].covariance)
    vnames = res[0].vparam_names
    cov = pd.DataFrame(_cov, columns=vnames, index=vnames)

    if not normalized:
        return cov

    for i, col in enumerate(cov.columns):
        cov[col] = cov[col]/stds[i]

    for i in range(len(cov)):
        cov.iloc[i] = cov.iloc[i]/stds[i]
    return cov
