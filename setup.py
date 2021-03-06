#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2008,2009,2014,2020 Jochen Küpper <jochen.kuepper@cfel.de>

import os
from setuptools import setup

extra_compile_args = []
library_dirs = []

long_description = """CMI Python extensions

Python extensions for calculations relevant to the manipulation of molecules.

Original author:    Jochen Küpper <jochen.kuepper@cfel.de>
Current maintainer: Jochen Küpper <jochen.kuepper@cfel.de>
"""


setup(name="cmiext",
      author              = "Jochen Küpper",
      author_email        = "jochen.kuepper@cfel.de",
      description         = "CMI Python extensions",
      license             = "GPL v3",
      url                 = "https://www.controlled-molecule-imaging.org/research/further_projects/software",
      version             = "1.0.1.dev0",
      long_description    = long_description,
      packages            = ['cmiext'],
      data_files          = [("", ["LICENSE", "COPYING", "README.devel.md"])]
      )
