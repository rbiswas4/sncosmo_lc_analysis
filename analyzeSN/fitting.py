#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import sncosmo
import simulate_lsst as sl
from numpy.random import normal, uniform


# Helper Functions to load the SNANA Data
def snanadatafile(snanafileroot, filetype='head', location='./'):
    '''
    obtain the name of the head of phot file of an SNANA simulation and dataset
    
    '''
    import os
    suffix = '_HEAD.FITS'
    if filetype == 'phot':
        suffix = '_PHOT.FITS'
    fname = snanafileroot + suffix
    return os.path.join(location, fname)
    
def loadSNANAData(snanafileroot, location='/.', snids=None, n=None):
    '''
    load a SNANA fits file into a list of `~astropy.Table.table` objects.
    
    
    Parameters
    ----------
    snanafileroot: string, mandatory
        root file name for the SNANA which is the prefix to '_HEAD.FITS', or '_PHOT.FITS'
    location: string, optional defaults to current working directory './' 
        directory where the head and phot files are located
    snids: integer/string, optional defaults to None
        if not None, only SN observations corresponding to SNID snid are loaded
    n: Integer, defaults to None
        if not None, only the first n SN light curves are loaded
        
        
    Returns: data
        list of `~astropy.Table.Table` each Table containing a light curve of a SN. 
        
    ..note: The column names of the SNANA data files are not reformated for SNCosmo use
    '''
    headfile = snanadatafile(snanafileroot, filetype='head', location=location)
    photfile = snanadatafile(snanafileroot, filetype='phot', location=location)
    data = sncosmo.read_snana_fits(head_file=headfile, phot_file=photfile, snids=snids, n=None)
    return data


# In[7]:

def addbands(sn, lsstbands, replacement):
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


def reformat_SNANASN(sn, lsstbands=None, replacements=None):
    '''
    reformat an SNANA light curve for use with SNCosmo
    
    Parameters
    ----------
    sn: `~astropy.Table.Table`, mandatory
        representing SNANA light curve
    lsstbands: list of strings, optional defaults to None
        list of unique strings in any of the 'FLT' column of SNANA files
    replacements: list of strings, optional defaults to None
        list of unique strings of the same size as lsstbands, and indexed in the 
        same order representing the keys in the sncosmo.bandpass registry for the
        same filters
    
    
    Returns
    -------
    `astropy.Table.Table` of the SNANA light curve reformatted for SNCosmo 
    '''
    #rename cols to names SNCosmo understands
    sn.rename_column("FLUXCAL",'flux')
    sn.rename_column("FLUXCALERR", 'fluxerr')
    #Add in SNANA magic ZP and sys
    sn["ZP"] = 27.5
    sn["ZPSYS"] = 'ab'
    # sn.rename_column('FLT', 'band')
    
    #Set up a truth dictionary from the metadata to set sim models
    truth ={}
    truth["c"] = sn.meta["SIM_SALT2c"]
    truth["x0"] = sn.meta["SIM_SALT2x0"]*10**(-0.4 * 0.27)
    truth["x1"] = sn.meta["SIM_SALT2x1"]
    truth["t0"] = sn.meta["SIM_PEAKMJD"]
    truth["mwebv"] = sn.meta["SIM_MWEBV"]
    truth["z"] = sn.meta["REDSHIFT_FINAL"]
    if replacements is not None:
        addbands(sn, lsstbands, replacements)
    return sn, truth



def lc(SNID, filename, location="./"):
    from astropy.io import fits
    summaryfile = location + '/' +filename + "_HEAD.FITS"
    photofile = location + '/' + filename + "_PHOT.FITS"

    
    summary = Table(fits.open(summaryfile)[1].data)
    
    snpointers = summary[np.array(map(int,summary["SNID"]))==SNID]
    ptrobs_min = snpointers["PTROBS_MIN"][0]
    ptrobs_max = snpointers["PTROBS_MAX"][0]
        
    lc = Table(fits.open(photofile)[1].data[ptrobs_min: ptrobs_max])
    lc.rename_column("FLUXCAL", "flux")
    lc.rename_column("FLUXCALERR","fluxerr")
    lc.rename_column("ZEROPT","ZP")
    #lc = photofile
    return lc

