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

"""Unit conversion routines"""

import const


def D2Cm(val):
    """Convert dipole moment from Debye to Coulomb * meter"""
    return val * 1e-21 / const.speed_of_light


def Cm2D(val):
    """Convert dipole moment from Coulomb * meter to Debye"""
    return val / (1e-21 / const.speed_of_light)


def Hz2J(val):
    """Hertz -> Joule"""
    return val * const.plancks_constant_h


def J2Hz(val):
    """Joule -> Hertz"""
    return val / const.plancks_constant_h


def V_m2kV_cm(val):
    """V/m -> kV/cm"""
    return val * 1e-5


def kV_cm2V_m(val):
    """kV/cm -> V/m"""
    return val / 1e-5
