# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008,2009 Jochen Küpper <software@jochen-kuepper.de>
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

"""Provide mathematical and physical constants.

Physical constants are based on CODATA values, see module codata.py for details.

Values last updated: $Date$"""

# mathematical constants
pi = 3.1415926535897931

# CODATA
from jkext.codata import codata
Boltzmann_constant                      = codata["Boltzmann constant"][0]
Planck_constant                         = codata["Planck constant"][0]
speed_of_light                          = codata["speed of light in vacuum"][0]
unified_atomic_mass                     = codata["unified atomic mass unit"][0]

# other physical units or conversion factors
Angstrom = 1e-10
