# -*- coding: utf-8 -*-
#
# Copyright (C) 2008 Jochen Küpper <software@jochen-kuepper.de>
# see LICENSE file for details

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

"""Provide physical constants.

If pyGSL is available, we will simply import its constants into our namespace.

Otherwise, we will provide "copies" of the most important constants,
using the same naming convention and the same values.

Values last updated: $Date$"""

try:
    from pygsl.const import *
except ImportError:
    angstrom = 1e-10
    pi = 3.1415926535897931
    plancks_constant_h = 6.6260689599999996e-34
    speed_of_light = 299792458.
    unified_atomic_mass = 1.6605387820000001e-27
    

### Local Variables:
### fill-column: 132
### End:
