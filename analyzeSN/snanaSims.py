#!/usr/bin/env python

"""
A module for helper functions to read SNANA simulations
"""

import sncosmo


class SnanaSims(object):
    """
    class to hold data from SNANA simulations and methods to manipulate the
    data


    Attributes
    ----------

    snList : list of `~astropy.Table.Table` with each Table containing a
        light curve of a SN. 
    
    """

    def __init__(self, headfile, photfile, snids=None, n=None):
        """
        Instantiate class from SNANA simulation output files in fits format

        Parameters
        ----------
        snids : integer/string, optional defaults to None
            if not None, only SN observations corresponding to SNID snid
            are loaded
        n : Integer, defaults to None
            if not None, only the first n SN light curves are loaded
        
        
        ..note: The column names of the SNANA data files are not reformated
                 for SNCosmo use
       """ 
        self.snList =  sncosmo.read_snana_fits(head_file=headfile,
                                               phot_file=photfile, 
                                               snids=snids, n=None)


    @classmethod
    def fromSNANAfileroot(cls, snanafileroot, location='./', snids=None,
        n=None):
        """
        Class constructor from a root file and a location

        Parameters
        ----------
        snanafileroot : string, mandatory
            root file name for the SNANA which is the prefix to
            '_HEAD.FITS', or '_PHOT.FITS'
        location : string, optional defaults to current working directory './' 
            directory where the head and phot files are located
        snids : integer/string, optional defaults to None
            if not None, only SN observations corresponding to SNID snid
            are loaded
        n : Integer, defaults to None
            if not None, only the first n SN light curves are loaded
        """
        headfile = self.snanadatafile(snanafileroot, filetype='head',
                                      location=location)
        photfile = self.snanadatafile(snanafileroot, filetype='phot',
                                      location=location)
        data = sncosmo.read_snana_fits(head_file=headfile,
                                       phot_file=photfile,
                                       snids=snids, n=None)
        return cls(head_file=headfile, phot_file=photfile, snids=snids,
                   n=n)


    @staticmenthod
    def snanadatafile(snanafileroot, filetype='head', location='./'):
        '''
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
            directory in which the file is located

        Returns
        -------
            string : absolute path to the SNANA file 

        '''
        import os
        suffix = '_HEAD.FITS'
        if filetype == 'phot':
            suffix = '_PHOT.FITS'
        fname = snanafileroot + suffix
        return os.path.join(location, fname)
 
    @staticmethod
    def addbandstoSN(sn, lsstbands, replacement):
        '''
        add a column called 'band' to the `~astropy.Table.Table` by 
        applying the map of lsstbands to replacements to the content
        of a column called 'FLT' 
        
        Parameters
        ----------
        sn: `~astropy.Table.Table` obtained by reading an SNANA light curve
        lsstbands: list of strings, mandatory
            list of strings representing the filters in sn, which can be found
            by `np.unique(sn['FLT'])
        replacements: list of strings, mandatory
            list of strings representing the filters as registered in SNCosmo in
            the same order as lsstbands
            
        Returns
        -------
        `~astropy.Table.Table` with 'FLT' column removed and 'band' column added
        '''
        filterarray = np.zeros(len(sn), dtype='S8')
        for i, flt in enumerate(lsstbands):
            mask = sn['FLT']==flt
            filterarray[mask] = replacement[i]
            band = Table.Column(filterarray, name='band', dtype='S8')
        sn.add_column(band)
        sn.remove_column('FLT')
    
    
