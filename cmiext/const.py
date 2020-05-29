# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of CMI Python extensions
# Copyright (C) 2008,2009,2010 Jochen Küpper <jochen.kuepper@cfel.de>
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
__author__ = "Jochen Küpper <jochen.kuepper@cfel.de>"

"""Provide mathematical and physical constants.

Physical constants are based on CODATA values as far as available, see module codata.py for details.
"""


# CODATA
from cmiext.codata import codata
atomic_unit_of_electric_dipole_moment   = codata["atomic unit of electric dipole mom."][0]
Bohr_radius                             = codata["Bohr radius"][0]
Boltzmann_constant                      = codata["Boltzmann constant"][0]
electron_charge                         = codata["elementary charge"][0]
electron_mass                           = codata["electron mass"][0]
electron_volt                           = codata["electron volt"][0]
Planck_constant                         = codata["Planck constant"][0]
speed_of_light                          = codata["speed of light in vacuum"][0]
unified_atomic_mass                     = codata["unified atomic mass unit"][0]
vacuum_permittivity                     = codata["electric constant"][0]
vacuum_permeability                     = codata["mag. constant"][0]


# other physical units or conversion factors
Angstrom                                = 1e-10 # m
Debye                                   = 1e-21 / speed_of_light # C m
inch                                    = 2.54e-2 # m
