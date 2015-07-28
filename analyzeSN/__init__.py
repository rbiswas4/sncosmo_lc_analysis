# from . import fitting
import os
from . import filters
from . import snanaSims

here = __file__
basedir = os.path.split(here)[0]
example_data = os.path.join(basedir, 'example_data')
