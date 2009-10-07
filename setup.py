#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2008,2009 Jochen Küpper <software@jochen-kuepper.de>


from numpy.distutils.core import setup, Extension
import sys

long_description = """JK Python extensions

Python extensions for calculations relevant to the manipulation of molecules.

Original authors:   Jochen Küpper <software@jochen-kuepper.de>
Current maintainer: Jochen Küpper <software@jochen-kuepper.de>
"""

version_major = sys.version_info[0]
version_minor = sys.version_info[1]

setup(name="jkext",
      author              = "Jochen Küpper",
      author_email        = "software@jochen-kuepper.de",
      description         = "JK Python extensions",
      license             = "GPL",
      url                 = "http://python.jochen-kuepper.de",
      version             = "0.2.1",
      long_description    = long_description,
      package_dir         = {'jkext': 'lib'},
      packages            = ['jkext'],
      ext_modules         = [Extension('jkext._wigner_gsl',
                                       sources = ['src/wigner_gsl.c', 'src/wigner_gsl_module.c'],
                                       libraries = ['gsl', 'gslcblas']),
                             Extension('jkext._wigner_avda',
                                       sources = ['src/wigner_avda.f',],
                                       extra_link_args = ['-bundle'],
                                       libraries = ['dl', 'python%d.%d' % (version_major,
                                                                           version_minor)]),
                             Extension('jkext._wigner_fft',
                                       sources = ['src/wigner_fft.f'],
                                       include_dirs = ['/opt/local/include'],
                                       extra_link_args = ['-bundle'],
                                       libraries = ['fftw3f', 'fftw3',
                                                    'dl', 'python%d.%d' % (version_major,
                                                                           version_minor)]),
                             ],
      scripts             = ['scripts/jkext_brute-force-orientation',
                             'scripts/jkext_calculate_energy',
                             'scripts/jkext_plot_energy',
                             'scripts/jkext_print_energy']
      )

