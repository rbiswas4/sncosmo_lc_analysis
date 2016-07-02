#!/usr/bin/env python

"""
A module for helper functions to read SNANA simulations
"""
from __future__ import absolute_import
import numpy as np
from astropy.io import fits
from astropy.utils import lazyproperty
from astropy.table import Table
import sncosmo
import pandas as pd

from . import filters

__all__ = ['SnanaSims']


class SnanaSims(object):

    """
    class to hold data from SNANA simulations and methods to manipulate the
    data


    Attributes
    ----------

    snList : list of `~astropy.table.Table`
        each Table contains a light curve of a SN. 

    """

    def __init__(self, headfile, photfile, snids=None, n=None, cutFunc=None):
        """
        Instantiate class from SNANA simulation output files in fits format

        Parameters
        ----------
        headfile : string, mandatory
            absolute path to the 'HEAD' file of the simulation
        photfile : string, mandatory
            absolute path to the 'HEAD' file of the simulation
        snids : integer/string, optional defaults to None
            if not None, only SN observations corresponding to SNID snid
            are loaded
        n : Integer, defaults to None
            if not None, only the first n SN light curves are loaded

        cutFunc : callable, optional, defaults to None
           function of the form cutFunc(snid, headfile, photfile) = True if the
           SN reperesented by snid in the files should be used. If None, no
           function is applied and consequently all SN in the simulation are
           used

        ..note: The column names of the SNANA data files are not reformated
                 for SNCosmo use


        """
        self.headfile = headfile
        self.headdata = self.get_headData(self.headfile)
        self.photfile = photfile
        self.phothdu = fits.open(photfile)

        self.snList = None # sncosmo.read_snana_fits(head_file=self.headfile,
        #  phot_file=self.photfile,
        #  snids=snids, n=None)

    @classmethod
    def fromSNANAfileroot(cls,
                          snanafileroot,
                          location='./',
                          snids=None,
                          n=None):
        """
        Class constructor from a root file and a location

        Parameters
        ----------
        snanafileroot : string, mandatory
            root file name for the SNANA which is the prefix to
            '_HEAD.FITS', or '_PHOT.FITS'
        location : string, optional defaults to current working directory './'
            Relative or absolute path to the directory where the head and phot
            files are located
        snids : integer/string, optional defaults to None
            if not None, only SN observations corresponding to SNID snid
            are loaded
        n : Integer, defaults to None
            if not None, only the first n SN light curves are loaded
        """

        headfile = cls.snanadatafile(snanafileroot, filetype='head',
                                     location=location)
        photfile = cls.snanadatafile(snanafileroot, filetype='phot',
                                     location=location)
        # data = sncosmo.read_snana_fits(head_file=headfile,
        #                               phot_file=photfile,
        #                               snids=snids, n=None)
        return cls(headfile=headfile, photfile=photfile, snids=snids,
                   n=n)

    @staticmethod
    def get_headData(headfile):
        """
        obtain the data in a headfile in a `pandas.DataFrame`

        Parameters
        ---------
        headfile : string, mandatory
            absolute path to an SNANA headFile

        Returns
        -------
        A dataFrame with the data in the Head file, 
        """
        headData = Table.read(headfile).to_pandas()
        # If the type is string (object in pandas)
        if headData.SNID.dtype == 'O':
            headData.SNID = headData.SNID.apply(lambda x: x.strip())
        return headData



    @staticmethod
    def snlc(snid, headData, phothdu):
        """
        obtain the SN light curve from the photometry files

        Parameters
        ----------
        headData : `panda.DataFrame`, mandatory
        phothdu : `hdu` 
        snid : index for SN. Can be integer or string.

        Returns
        -------
        `numpy.recarray` of data
        """
        if not isinstance(snid, int):
            # snid is a string
            snid = snid.strip().lower()
        
        ptrs_cols = ['PTROBS_MIN', 'PTROBS_MAX']
        ptrs = headData.query('SNID == @snid')[ptrs_cols].values
        ptrs[0] -= 1

        data = phothdu[1].data
        dt = data.dtype

        lc = []
        for i, datum in enumerate(data):
            if i >= ptrs[0]:
                lc.append(datum)
                if i > ptrs[1]:
                    break
        return np.array(lc, dtype=dt)


    @staticmethod
    def snanadatafile(snanafileroot, filetype='head', location='./'):
        """
        obtain the name of the head or phot file of an SNANA simulation
        and dataset

        Parameters
        ----------
        snanafileroot : string, mandatory
            root file name for the SNANA which is the prefix to
            '_HEAD.FITS', or '_PHOT.FITS'
        filetype : string, optional defaults to 'head'
            'head' or 'phot' depending on whether a summary file or a photometry
            file is being used.
        location : string, optional defaults to current working directory './' 
            relative or absolute path to the directory in which the file is
            located

        Returns
        -------
            string : absolute path to the SNANA file 

        """
        import os

        desiredfiletype = ['head', 'phot']
        filetype = filetype.lower()
        if not filetype in desiredfiletype:
            raise ValueError(
                'filetype should be one of "head" or "phot"', filetype)
        location = os.path.abspath(location)
        suffix = '_HEAD.FITS'
        if filetype.lower() == 'phot':
            suffix = '_PHOT.FITS'
        fname = snanafileroot + suffix
        return os.path.join(location, fname)

    @staticmethod
    def addbandstoSN(sn, snanaBands, replacement):
        '''
        add a column called 'band' to the `~astropy.Table.Table` by 
        applying the map of lsstbands to replacements to the content
        of a column called 'FLT' 

        Parameters
        ----------
        sn: `~astropy.Table.Table` obtained by reading an SNANA light curve
        snanaBands: list of strings, mandatory
            list of strings representing the filters in sn, which can be found
            by `np.unique(sn['FLT'])
        replacements: list of strings, mandatory
            list of strings representing the filters as registered in SNCosmo in
            the same order as lsstbands

        Returns
        -------
        `~astropy.Table.Table` with 'FLT' column removed and 'band' column added
        '''

        from astropy.table import Table

        filterarray = np.zeros(len(sn), dtype='S8')
        for i, flt in enumerate(snanaBands):
            mask = sn['FLT'] == flt
            filterarray[mask] = replacement[i]
            band = Table.Column(filterarray, name='band', dtype='S8')
        sn.add_column(band)
        sn.remove_column('FLT')

    @staticmethod
    def reformat_SNANASN(sn, snanaBands=None, replacements=None):
        '''
        reformat an SNANA light curve for use with SNCosmo

        Parameters
        ----------
        sn: `~astropy.Table.Table`, mandatory
            representing SNANA light curve
        snanaBands: list of strings, optional defaults to None
            list of unique strings in any of the 'FLT' column of SNANA files
        replacements: list of strings, optional defaults to None
            list of unique strings of the same size as lsstbands, and indexed in the 
            same order representing the keys in the sncosmo.bandpass registry for the
            same filters


        Returns
        -------
        `astropy.Table.Table` of the SNANA light curve reformatted for SNCosmo 
        '''

        from astropy.table import Table

        # rename cols to names SNCosmo understands
        sn.rename_column("FLUXCAL", 'flux')
        sn.rename_column("FLUXCALERR", 'fluxerr')
        # Add in SNANA magic ZP and sys
        sn["ZP"] = 27.5
        sn["ZPSYS"] = 'ab'

        if replacements is not None:
            SnanaSims.addbandstoSN(sn, snanaBands, replacements)
        else:
            sn.rename_column('FLT', 'band')
        return sn

    @staticmethod
    def matchSNANAbandnamesinregistry():
        """
        Will have to build this along as we go, as I don't know the variety
        of naming conventions

        """

        bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
        for bandpass in bandPassList:
            band = sncosmo.get_bandpass('LSST_' + bandpass)
            band.name = bandpass
            if bandpass == 'y':
                band.name = 'Y'
            sncosmo.registry.register(band)
