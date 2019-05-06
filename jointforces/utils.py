import os
import sys
import numpy as np
from distutils.sysconfig import get_python_lib


def set_gmsh_path(path, overwrite=False):
    # find site-package directory
    dir = get_python_lib()
    pathfile = dir + '/gmsh.pth'

    # add lib subdir
    path = path + '/lib'

    if os.path.exists(pathfile) and not overwrite:
        print('Gmsh path already set.')
    else:
        # write path file
        with open(pathfile, 'w+') as f:
            f.write(path+'\n')

        # add it now to eliminate the need to restart Python
        sys.path.append(path)
        print('Set Gmsh path to: {}'.format(path))


def set_saeno_path(path, overwrite=False):
    # find site-package directory
    dir = get_python_lib()
    pathfile = dir + '/saeno.cfg'

    if os.path.exists(pathfile) and not overwrite:
        print('SAENO path already set.')
    else:
        # write path file
        with open(pathfile, 'w+') as f:
            f.write(path + '\n')

        # add it now to eliminate the need to restart Python
        os.environ["PATH"] += os.pathsep + path

        # add it now to eliminate the need to restart Python
        print('Set SAENO path to: {}'.format(path))


def load(file):
    return np.load(file).item()
