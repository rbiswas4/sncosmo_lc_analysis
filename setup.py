#from ez_setup import use_setuptools
#use_setuptools()
from setuptools import setup, find_packages

setup(# package information
      name="analyzeSN",
      version="0.0.1dev",
      description='A set of utilities to analyze SN light curves',
      long_description=''' ''',
      # What code to include as packages
      packages=['analyzeSN'],
      # What data to include as packages
      include_package_data=True,
      package_data={'': ['example_data/*.FITS']}

      )
