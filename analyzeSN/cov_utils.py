#!/usr/bin/env python
"""
A number of utility functions to conventientyly deal with covariances
"""

import numpy as np
import pandas as pd
# from copy import deepcopy


def log_covariance(covariance, paramName, paramValue, factor=1.):
    """
    Covariance of the parameters with parameter paramName replaced by
    factor * np.log(param) everywhere, and its true value is paramValue,
    assuming linear propagation

    Parameters
    ----------
    covariance : `pandas.DataFrame`, mandatory
        representing covariance matrix
    paramName : int or str, mandatory
        integer or parameter name specifying the position of the variable
        whose logarithm must be taken
    paramValue : float, mandatory
        true/estimated value of the variable itself
    factor : float, optional, defaults to 1.
        Factor multiplying the natural logarithm. For example,
        if the relevant transformation is going from 'f' to
        -2.5 log10(f), the factor should be -2.5 /np.log(10)

    Returns
    -------

    Examples
    --------
    """
    if isinstance(paramName, np.int):
        cov = covariance.values
        cov[:, paramName] = factor * cov[:, paramName] / paramValue
        cov[paramName, :] = factor * cov[paramName, :] / paramValue
        return cov

    covariance[paramName] = factor * covariance[paramName] / paramValue
    covariance.loc[paramName] = factor * covariance.loc[paramName] / paramValue

    return covariance


def subcovariance(covariance, paramList, array=False):
    """
    returns the covariance of a subset of parameters in a covariance dataFrame.

    Parameters
    ----------
    covariance : `pandas.DataFrame` representing square covariance matrix
        with parameters as column names, and index as returned by covariance
    paramList : list of strings, mandatory
        list of parameters for which the subCovariance matrix is desired.
        The set of parameters in paramList must be a subset of the columns
        and indices of covariance
    array : boolean, optional, defaults to False
        if true, return `numpy.ndarray`, if False return `pandas.DataFrame`
    Returns
    -------
    """
    df = covariance.ix[paramList, paramList]
    if array:
        return df.values
    else:
        return df


def covariance(covArray, paramNames=None, normalized=False):
    """
    Returns the covariance Matrix for the varied parameters as a
    `pandas.DataFrame`, so that covariance elements may be called
    by either index or parameters.

    Parameters
    ----------
    covArray : `numpy.ndarray` of the covariance, mandatory
    paramNames : list of strings, optional, defaults to None
    normalized : Bool, optional, defaults to False
        whether to return the normalized covariance matrix

    Returns
    -------
    a `pandas.DataFrame` with column names and indexes given by the parameter
    names. If paramNames is None, the return is a DataFrame with indexes and
    column names chosen by pandas.

    Examples
    --------
    >>> cov = covariance(covArray, paramNames=['t0', 'x0', 'x1, 'c'])
    >>> cov.ix[['t0', 'x1'],['t0', 'x1']]
    >>> cov.iloc[[0, 2], [0, 2]]
    """

    l, w = np.shape(covArray)
    # Check for the covariance matrix being square, not checking for symmetry
    if l != w:
        raise ValueError('The covariance matrix is not square; length!=width')

    if paramNames is not None:
        if len(paramNames) != w:
            raise ValueError('The number of parameters must match the length'
                             ' of the covariance matrix')
        cov = pd.DataFrame(covArray, columns=paramNames, index=paramNames)
    else:
        cov = pd.DataFrame(covArray)

    if not normalized:
        return cov

    # normalize if requested
    stds = cov.values.diagonal()
    for i, col in enumerate(cov.columns):
        cov[col] = cov[col]/stds[i]

    for i in range(len(cov)):
        cov.iloc[i] = cov.iloc[i]/stds[i]
    return cov
