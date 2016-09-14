# IPython log file

import analyzeSN as ans
import os
def test_load():
    headFile = os.path.join(ans.__path__[0], 'example_data', 'snana_fits_HEAD.FITS')
    photFile = os.path.join(ans.__path__[0], 'example_data', 'snana_fits_PHOT.FITS')
    sne = ans.SNANASims(headFile=headFile, photFile=photFile, coerce_inds2int=False)
    assert len(sne.headData) == 2
    assert len(sne.get_SNANA_photometry(snid='03d1aw')) > 0
    assert len(sne.get_SNANA_photometry(snid='03d1ax')) > 0
def test_snanadatafiles():
    headFile = os.path.join(ans.__path__[0], 'example_data', 'snana_fits_HEAD.FITS')
    photFile = os.path.join(ans.__path__[0], 'example_data', 'snana_fits_PHOT.FITS')
    #sne = ans.SNANASims(headFile=headFile, photFile=photFile, coerce_inds2int=False)
    loc = os.path.join(ans.__path__[0], 'example_data')
    testheadFile = ans.SNANASims.snanadatafile(snanafileroot='snana_fits',
                                               filetype='head',
                                               location=loc)
    testphotFile = ans.SNANASims.snanadatafile(snanafileroot='snana_fits',
                                               filetype='pHot',
                                               location=loc)
    assert testheadFile == headFile
    assert testphotFile == photFile

