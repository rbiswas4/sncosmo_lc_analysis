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
