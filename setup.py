#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) Gabriele Vajente (2016)

"""Setup the BruCo package
"""

import glob
import os.path

try:
    import setuptools
except ImportError:
    import ez_setup
    ez_setup.use_setuptools()
finally:
    from setuptools import (setup, find_packages)

# set basic metadata
PACKAGENAME = 'bruco'
DISTNAME = 'bruco'
AUTHOR = 'Duncan Macleod'
AUTHOR_EMAIL = 'duncan.macleod@ligo.org'
LICENSE = 'GPLv3'

# Use the find_packages tool to locate all packages and modules
packagenames = find_packages()

# glob for all scripts
if os.path.isdir('bin'):
    scripts = glob.glob(os.path.join('bin', '*'))
else:
    scripts = []

# -----------------------------------------------------------------------------
# run setup

setup(name=DISTNAME,
      provides=[PACKAGENAME],
      version=None,
      description=None,
      long_description=None,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      packages=packagenames,
      include_package_data=True,
      scripts=scripts,
      requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'glue',
          'pylal',
      ],
      install_requires=[
      ],
      dependency_links=[
          'http://software.ligo.org/lscsoft/source/glue-1.49.1.tar.gz',
      ],
      use_2to3=True,
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics',
          'Operating System :: POSIX',
          'Operating System :: Unix',
          'Operating System :: MacOS',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      ],
)
