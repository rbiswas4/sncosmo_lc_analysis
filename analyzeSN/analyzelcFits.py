#!/usr/bin/env python
"""
classes and methods to analyze the outputs of light curve characterization
routines, particularly for SALT2 light curves
"""

import numpy as np
import pandas as pd
# from copy import deepcopy
from . import cov_utils as cutils

class ResChar(object):
    """
    A class to hold results of characterizing a light curve fit. Mostly
    this will be constructed from `sncosmo.utils.res` instances
    """

    def __init__(self,
                 vparam_names,
                 param_names,
                 parameters,
                 covariance,
                 samples=None,
                 weights=None,
                 sncosmoModel=None):
        """
        constructor for class

        Parameters
        ----------
        vparam_names : list of strings, mandatory
            model parameters inferred in the characterization
        param_names : list of strings, mandatory
            model parameters (complete list)
        parameters : list of floats, mandatory
            values of all model parameters in same order as param_names
        covariance : `numpy.ndarray`, mandatory
            Covariance of all the varied parameters. Should be a
            symmetric positive definite array of shape(n, n) where
            n = len(vparam_names)
        samples : `numpy.ndarray`, optional, defaults to None
            samples of shape(num, n) where num = number of samples, and
            n = len(vparam_names). May not be independent
        weights : `np.array` , optional defaults to None
            must have length equal to len(samples) if present. For
            mcmc_lc, this is usually used as `np.ones`. For nest_lc
            the weights are used to calculate quantities
        sncsomoModel : `sncosmo.Model`, optional, defaults to None
            model returned from sncsomo estimate
        """

        self.vparam_names = vparam_names
        self.param_names = param_names
        self._parameters = parameters
        self._covariance = covariance
        self.samples = samples
        self.weights = weights
        self.sncosmoModel = sncosmoModel


    @classmethod
    def fromSNCosmoRes(cls, SNCosmoRes):
        """
        Instantiate this object from an instance of `sncosmo.utils.res`

        Parameters
        ----------
        SNCosmoRes : instance of `sncosmo.utils.res
        """

        res, model = SNCosmoRes

        # samples if the method was mcmc/nest_lc but not if max_lc
        # weights makes sense for mcmc methods
        samples = None
        if 'samples' in res.keys():
            samples = res['samples']
            weights = np.ones(len(samples))

        # if method was nest_lc
        if 'weights' in res.keys():
            weights = res['weights']

        return cls(vparam_names=res.vparam_names,
                   param_names=res.param_names,
                   parameters=res.parameters,
                   covariance=res.covariance,
                   samples=samples,
                   weights=weights,
                   sncosmoModel=model)



    @property
    def parameters(self):
        """
        return the model parameters as a `pd.Series`
        """
        return pd.Series(self._parameters, index=self.param_names)
    @property
    def vparams(self):
        """
        return the values of the varied parameters as a `pd.Series`
        """
        vparameters = [self._parameters[self.param_names.index(v)]
                       for v in self.vparam_names]
        return pd.Series(vparameters, index=self.vparam_names)

    @property
    def covariance(self):
        """
        return the covariance as a `pd.DataFrame`
        """
        cutils.covariance(self._covariance, paramNames=self.vparam_names)

