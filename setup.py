#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2008 Jochen Küpper <software@jochen-kuepper.de>


import os
from distutils.core import setup

extra_compile_args = []

# library dirs
library_dirs = []


### No changes needed below ###

long_description = """JK Python extensions

Contains many useful extension I came across over the time.

Original authors:   Jochen Küpper <software@jochen-kuepper.de>
Current maintainer: Jochen Küpper <software@jochen-kuepper.de>
"""


setup(name="beamline",
      author              = "Jochen Küpper",
      author_email        = "software@jochen-kuepper.de",
      description         = "JK Python extensions",
      license             = "GPL",
      url                 = "http://python.jochen-kuepper.de",
      version             = "0.0.1",
      long_description    = long_description,      
      package_dir         = {'jkext': 'lib'},
      packages            = ['jkext'],
      )

