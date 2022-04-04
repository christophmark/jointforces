#!/usr/bin/env python

from setuptools import setup

setup(
    name='jointforces',
    packages=['jointforces'],
    version='1.0.0',
    description='A Python package to conduct 3D Traction Force Microcopy on multicelluar aggregates (spheroids).',
    url='',
    download_url = '',
    author='Christoph Mark, David Böhringer',
    author_email='christoph.mark@fau.de',
    license='The MIT License (MIT)',
    install_requires=['saenopy>=0.7.1', 
                      'gmsh-sdk>=4.5.0',
                      'numpy>=1.16.2',
                      'pandas>=1.2.0',
                      'matplotlib>=2.2.2',
                      'scipy>=1.1.0',
                      'scikit-image>=0.14.2',
                      'openpiv>=0.23.4',
                      'tqdm>=4.26.0',
                      'natsort>=5.1.0', 
                      'dill>=0.2.9',
                      'roipoly>=0.5.2',
                      'openpyxl',
                      'pyyml',
					  'matplotlib-scalebar'],
    keywords = ['piv', 'contractility', 'material simulation', 'biophysics'],
    classifiers = [],
    )
