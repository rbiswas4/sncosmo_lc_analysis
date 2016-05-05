from distutils.core import setup

setup(# package information
      name="analyzeSN",
      version="0.0.2dev",
      description='A set of utilities to analyze SN light curves',
      long_description=''' ''',
      # What code to include as packages
      packages=['analyzeSN'],
      package_dir={'analyzeSN':'analyzeSN'},
      # What data to include as packages
      include_package_data=True,
      package_data={'analyzeSN': ['example_data/*.FITS', 'example_data/*.dat', 'example_data/*.md']}
      )
