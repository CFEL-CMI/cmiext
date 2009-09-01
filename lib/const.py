# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008 Jochen Küpper <software@jochen-kuepper.de>
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this programm for scietific work, you must correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.
from __future__ import division

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
    boltzmann = 1.3806504000000001e-23
    pi = 3.1415926535897931
    plancks_constant_h = 6.6260689599999996e-34
    speed_of_light = 299792458.
    unified_atomic_mass = 1.6605387820000001e-27
    vacuum_permittivity = 8.854187817e-12
