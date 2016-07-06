#!/usr/bin/env python

"""
A module for helper functions to read SNANA simulations
"""
from __future__ import absolute_import
import numpy as np
from astropy.io import fits
import fitsio
from astropy.utils import lazyproperty
from astropy.table import Table, Column
import sncosmo
import pandas as pd
from collections import OrderedDict as odict
from . import filters

__all__ = ['SnanaSims']


class SnanaSims(object):

    """
    class to hold data from SNANA simulations and methods to manipulate the
    data.


    Attributes
    ----------

    snList : list of `~astropy.table.Table`
        each Table contains a light curve of a SN. 

    """

    def __init__(self, headfile, photfile, snids=None, n=None, cutFunc=None,
                 noPhot=False):
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
           function of the form cutFunc(snid, headDict, photData=None) which
           should return True if a SN with SNID snid, headData represented by
           a `collections.orderedDict` instance headDict, and photData
           represented by a `numpy.recarray` are passed
           SN reperesented by snid in the files should be used. If None, no
           function is applied and consequently all SN in the simulation are
           used. If noPhot is set to True,then the cutFunc can be applied on
           the headData only
        noPhot : Bool, optional, defaults to False
            If True, it is assumed that the cuts are only applied on the
            headData

        ..note: The column names of the SNANA data files are not reformated
                 for SNCosmo use


        """
        self.headfile = headfile
        self.headData = self.get_headData(self.headfile)
        self.photfile = photfile
        self.phothdu = fits.open(photfile)
        self.photFits = fitsio.FITS(photfile)
        self.snids = None
        self._snList = None # sncosmo.read_snana_fits(head_file=self.headfile,

        @property
        def snList(self):
            """
            Return a list of light curves of selected SN with the selection
            determined by self.snids. 

            If self.snids is None, all SN in the simulation shall be accessed.
            If one only wants to access a subset of the SN, the SnanaSims object
            should be instantiated, and a suitable selection performed and
            assigned to self.snids before accessing this attribute
            """
            if self._snList is None:
                return self.get_snList(self.snids)
            elif len(set(self.snids) - set(self.snidsForsnList(self._snList)))\
                > 0 :
                newids = \
                list(set(self.snids) - set(self.snidsForsnList(self._snList)))
                newSN = self.get_snList(newids)
                return self._snList + self.newSN
            else:
                return self._snList

    @staticmethod
    def snidsForSnList(snlist):
        """
        return a list of snids for each SN in SNList
        """
        return list(sn.meta['SNID'] for sn in snlist)


#    def getSN(self, cutFunc, noPhot=False, numstart=0, nums=10000):
#
#        NoSelections = False
#        select = False
#        if cutFunc is None:
#            NoSelections = True
#        if nums is None:
#            headData = self.headData.reset_index()
#        else:
#            headData = self.headData.reset_index().head(nums+numstart)
#
#	for ind, row in headData.iterrows():
#            if ind < numstart:
#                continue
#	    summaryProps = odict(row)
#	    snid = summaryProps['SNID']
#	    if noPhot:
#		photData = None
#	    else:
#		photData = self.get_SNANAPhotometry(snid, self.headData,
#						    self.phothdu)
#            if not NoSelections:
#                select = cutFunc(snid, summaryProperties=summaryProps,
#                                  photometryData=photData, noPhot=True)
#            
#            if select or NoSelections:
#		if photData is None:
#		    lc = sne.get_SNANAPhotometry(snid, sne.headData,
#					         sne.phothdu)
#		else:
#		    lc = photData
#		myTable = Table(lc)
#		myTable.meta = summaryProps
#		yield myTable
    @staticmethod
    def addColumnsToLC(snid, format='pandas'):

        lc = self.get_SNANAPhotometry(snid,
                                      headData=self.headData
                                      photFITS=self.photFits)

        lct = pd.DataFrame(lc)
        lct['SNID'] = snid
        return lct

    def getAllSNPhot(self, cutFunc=None, noPhot=True, snids=None, range=(0, -1)):
        """
        """
        if snids is None
            snids = self.headData.index.values[range[0]:range[1]]
        return (self.addColumnsToLC(snid) for snid in snids if cutFunc(snid,
                                                                       headData,
                                                                       self.addColumnsToLC))

    def getSelectedSNPhot(self, cutFunc=None, noPhot=True, snids=None, range=(0, -1)):
        """
        """
        if snids is None
            snids = self.headData.index.values[range[0]:range[1]]
        return (self.addColumnsToLC(snid) for snid in snids)

    def getSN(self, cutFunc=None, noPhot=True, photPandas=False):
        """
        cutFunc
        noPhot :
        photPandas:
        """
        headData = self.headData.reset_index()
	for ind, row in headData.iterrows():
	    summaryProps = odict(row)
	    snid = summaryProps['SNID']
	    lc = self.get_SNANAPhotometry(snid,
                                          headData=self.headData,
	          		          photFits=self.photFits)
	    myTable = Table(lc)
	    myTable.meta = summaryProps
	    yield myTable
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
        headData = fitsio.read(headfile)
        # If the type is string
        if headData['SNID'].dtype.type is np.string_:
            data = headData['SNID'].data
            name = headData['SNID'].name
            dtype = headData['SNID'].dtype
            arr = list(x.strip().lower() for x in data)
            col = Column(data=arr, name=name, dtype=dtype) 
            headData['SNID'] = col
        headData = headData.to_pandas()
        return headData.set_index('SNID')



    @staticmethod
    def get_SNANAPhotometry(snid,
                            headData=None,
                            photFits=None,
                            headFile=None,
                            photFile=None):
        """
        obtain the SNANA photmetry from the head and photometry
        files

        Parameters
        ----------
        snid : string or int, mandatory
            index of SNANA SN, If string, whitespace padding will be removed
            and all characters are transformed into lower case.
        headData : `panda.DataFrame`, optional, defaults to None
            a `pandas.DataFrame` object holding data from the SNANA
            headfile. If None, uses the headFile paramter to generate this
            and headFile cannot be None in that case.
        photFits : `fitsio.FITS` instance, optional, defaults to None
            `fitsio.FITS` instance corresponding to the SNANA photFile, If
            None, generates this from the `photFile` parameter
        #    from the photFile parameter and photFile cannot be None.
        headFile : string, optional, defaults to None
            absolute path to SNANA photFile. Ignored if headData is supplied
        photFile : string, optional, defaults to None
            absolute path to SNANA photFile. Ignored if phothdu is supplied

        Returns
        -------
        `numpy.recarray` of data
        """
        if isinstance(snid, basestring):
            # snid is a string
            snid = snid.strip().lower()
        
        if headData is None:
            if headFile is None:
                raise ValueError('both headData and headFile cannot be None\n')
            headData = SnanaSims.get_headData(headFile) 

        ptrs_cols = ['PTROBS_MIN', 'PTROBS_MAX']
        ptrs = headData.ix[snid, ptrs_cols].values
        ptrs[0] -= 1

        # if the `fitsio.FITS` instance is not provided
        if photFits is None:
            if photFile is None:
                raise ValueError('both headData and headFile cannot be None\n')
            photFits = fitsio.FITS(photFile)

        lcdata = photFits[1][ptrs[0]: ptrs[1]]
        return lcdata



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
            string : absolute path to the SNANA files 

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
