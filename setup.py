#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, unicode_literals

import os
import re
from io import open  # Python 2
from os import path

from setuptools import find_packages, setup

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


def lookup_version():
    with open(os.path.join('zernike', 'version.py'), 'r') as f:
        m = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", f.read(), re.M)
    return m.group(1)


setup(
    name='zernike',
    version=lookup_version(),
    description='Python code for Zernike polynomials',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/jacopoantonello/zernike',
    author='Jacopo Antonello',
    author_email='jacopo@antonello.org',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3',
    ],
    install_requires=['numpy', 'h5py'],
    extras_require={'plotting': ['matplotlib']},
    packages=find_packages(exclude=['tests*', 'examples*']),
    python_requires='>=2.7',
)
