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
import numpy


def D2Cm(val):
    """Convert dipole moment from Debye to Coulomb * meter"""
    return numpy.array(val) * 1e-21 / const.speed_of_light


def Cm2D(val):
    """Convert dipole moment from Coulomb * meter to Debye"""
    return numpy.array(val) / (1e-21 / const.speed_of_light)


def Hz2J(val):
    """Hertz -> Joule"""
    return numpy.array(val) * const.plancks_constant_h


def MHz2J(val):
    """Mega-Hertz -> Joule"""
    return Hz2J(val * 1e6)


def invcm2J(val):
    """cm^{-1} -> Joule"""
    return val * 100 * const.plancks_constant_h * const.speed_of_light


def J2Hz(val):
    """Joule -> Hertz"""
    return numpy.array(val) / const.plancks_constant_h


def J2MHz(val):
    """Joule -> Mega-Hertz"""
    return J2Hz(val) / 1e6


def J2invcm(val):
    """Joule -> cm^{-1}"""
    return val / const.plancks_constant_h / const.speed_of_light / 100


def V_m2kV_cm(val):
    """V/m -> kV/cm"""
    return numpy.array(val) * 1e-5


def kV_cm2V_m(val):
    """kV/cm -> V/m"""
    return numpy.array(val) / 1e-5

def A32CM2_V(val):
    """Å -> CM^2/V"""
    return val*const.vacuum_permittivity*((const.angstrom)**3)*4*const.pi

def dcfields2omega(dcfields,rotcon,dipole):
    
    return numpy.array(dcfields)*dipole/rotcon

def acfields2deltaomega(acfields,rotcon,polarizability):

    return acfields**2*polarizability/(4*rotcon)

def omega2dcfields(omega,rotcon,dipole):
    assert dipole != 0
    return numpy.array(omega)*rotcon/dipole

def deltaomega2acfields(deltaomega,rotcon,polarizability):
    assert polarizability != 0
    return (numpy.array(deltaomega)*4*rotcon/polarizability)**0.5
