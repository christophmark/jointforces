# check if SAENO path is stored
import os
from distutils.sysconfig import get_python_lib

cfg_file = get_python_lib() + '/saeno.cfg'

if os.path.exists(cfg_file):
    with open(cfg_file, 'r') as f:
        SAENOPATH = f.readlines()[0].strip()
else:
    SAENOPATH = '.'

# import package files
from . import mesh
from . import materials
from . import simulation
from . import piv
from . import force
from . import experiment
from .utils import set_gmsh_path, set_saeno_path, load
