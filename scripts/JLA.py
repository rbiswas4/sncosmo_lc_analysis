import numpy as np
from astropy.table import Table


class PhotometryData(object):
    def __init__(self, fname):
        self.fname = fname

        self.band_prefix = 'SNLS_'
        with open(self.fname) as f:
            contents = f.read()
        headcontents, photcontents = contents.split('#end\n')
        self.contents = contents
        self.headcontents = headcontents
        self.photcontents = photcontents
        self.alias = {'time': 'date', 'band': 'filter', 'zpsys':'magsys'}

    @property
    def names(self):
        _headers =  self.headers()

        names = []
        for col in _headers:
            tmp = col.lower()
            for key, vals in self.alias.items():
                if col.lower() in vals:
                    tmp = key
            names.append(tmp.lower())
        return names


    def headers(self):
        names = [x.strip()[1:-1] for x in self.headcontents.split('\n')
                 if x.startswith('#')]
        return names

    def photometryL(self):
        photList = [x.strip().split() for x in self.photcontents.split('\n')]
        photTransposed = [x for x in photList if x !=[]]
        return [x for x in zip(*photTransposed)]

    @property
    def meta(self):
        lines = [x.strip()[1:].split() for x in self.headcontents.split('\n')
                 if x.startswith('@')]
        keys, values = zip(*lines)
        mydict = {}
        for pairs in lines:
            mydict[pairs[0]] = pairs[1]
        return mydict

    @staticmethod
    def sanitize(table, bandprefix):
        table['zpsys'] = 'ab'
        table['band'] = [bandprefix + x[-1] for x in np.asarray(table['band'])]
        return table

        
    @property
    def photometryTable(self):
        pT =  Table(self.photometryL(), names=self.names)
        pT = self.sanitize(pT, self.band_prefix)
        for var in ['time', 'flux', 'fluxerr', 'zp']:
            _s = pT[var].astype(np.float) 
            pT.replace_column(var, _s)
        pT.meta = self.meta
        return pT


