#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, unicode_literals

import os
import re
import sys
from io import open  # Python 2
from os import path
from subprocess import check_output

from setuptools import find_packages, setup

# Get the long description from the relevant file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


def update_version():
    try:
        toks = check_output('git describe --tags --long --dirty',
                            universal_newlines=True,
                            shell=True).strip().split('-')
        version = toks[0].strip('v') + '.' + toks[1]
        last = check_output('git log -n 1',
                            universal_newlines=True,
                            shell=True)
        date = re.search(r'^Date:\s+([^\s].*)$', last, re.MULTILINE).group(1)
        commit = re.search(r'^commit\s+([^\s]{40})', last,
                           re.MULTILINE).group(1)
    except Exception as e:
        version = '0.0.0'
        date = 'unknown'
        commit = 'unknown'
        print('Cannot update version {}'.format(str(e)), file=sys.stderr)

    with open(path.join('zernike', 'version.py'), 'w', newline='\n') as f:
        f.write('#!/usr/bin/env python\n')
        f.write('# -*- coding: utf-8 -*-\n\n')
        f.write("__version__ = '{}'\n".format(version))
        f.write("__date__ = '{}'\n".format(date))
        f.write("__commit__ = '{}'".format(commit))


update_version()


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
    license='GPLv3+',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Physics',
        ('License :: OSI Approved :: GNU General Public License v3 ' +
         'or later (GPLv3+)'),
        'Programming Language :: Python :: 3',
    ],
    setup_requires=['numpy', 'h5py'],
    extras_requires={'plotting': ['matplotlib']},
    packages=find_packages(exclude=['tests*', 'examples*']),
    python_requires='>=2.7',
)
