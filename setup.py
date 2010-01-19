#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2008,2009 Jochen Küpper <software@jochen-kuepper.de>


import os
from distutils.core import setup

extra_compile_args = []
library_dirs = []

long_description = """JK Python extensions

Python extensions for calculations relevant to the manipulation of molecules.

Original authors:   Jochen Küpper <software@jochen-kuepper.de>
Current maintainer: Jochen Küpper <software@jochen-kuepper.de>
"""


setup(name="jkext",
      author              = "Jochen Küpper",
      author_email        = "software@jochen-kuepper.de",
      description         = "JK Python extensions",
      license             = "GPL",
      url                 = "http://python.jochen-kuepper.de",
      version             = "0.3.0",
      long_description    = long_description,
      package_dir         = {'jkext': 'lib'},
      packages            = ['jkext'],
      scripts             = ['scripts/jkext_GAMESS',
                             'scripts/jkext_linearize']
      )

