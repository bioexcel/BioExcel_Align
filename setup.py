#! /usr/bin/env python

'''
Setup script for bioexcel_align.
'''

from distutils.command.install import INSTALL_SCHEMES
from setuptools import setup

for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']

setup(
    name='bioexcel_align',
    version='0.2.0',
    description=('Alignment workflow python package'),
    author='Darren White',
    author_email='d.white@epcc.ed.ac.uk',
    scripts=['bin/bxcl_align'],
    packages=['bioexcel_align'],
    package_dir={'bioexcel_align': 'bioexcel_align'}
)
