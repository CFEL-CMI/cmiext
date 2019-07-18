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


__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

"""Unit conversion routines"""

import cmiext.const
import cmiext.codata
import numpy


def D2Cm(val):
    """Convert dipole moment from Debye to Coulomb * meter"""
    return numpy.array(val) * 1e-21 / cmiext.const.speed_of_light


def Cm2D(val):
    """Convert dipole moment from Coulomb * meter to Debye"""
    return numpy.array(val) / (1e-21 / cmiext.const.speed_of_light)


def eV2m(val):
    """eV -> wavelenght (m)"""
    return 1 / numpy.array(val) / cmiext.codata.codata["electron volt-inverse meter relationship"][0]


def E2I(val):
    """field amplitude (V/m) -> intensity (m)"""
    return 0.5 * cmiext.const.speed_of_light * cmiext.const.vacuum_permittivity * numpy.array(val)**2


def Hz2J(val):
    """Hertz -> Joule"""
    return numpy.array(val) * cmiext.const.Planck_constant


def MHz2J(val):
    """Mega-Hertz -> Joule"""
    return Hz2J(val * 1e6)


def I2E(val):
    """intensity (W/m^2) -> amplitude (V/m)"""
    return numpy.sqrt(2 * val / (cmiext.const.speed_of_light * cmiext.const.vacuum_permittivity))


def invcm2J(val):
    """cm^{-1} -> Joule"""
    return val * 100 * cmiext.const.Planck_constant * cmiext.const.speed_of_light


def J2eV(val):
    """Joule -> electron volt"""
    return cmiext.convert.numpy.array(val) / cmiext.codata.codata["electron volt-joule relationship"][0]


def J2Hz(val):
    """Joule -> Hertz"""
    return numpy.array(val) / cmiext.const.Planck_constant


def J2MHz(val):
    """Joule -> Mega-Hertz"""
    return J2Hz(val) / 1e6


def J2invcm(val):
    """Joule -> cm^{-1}"""
    return val / cmiext.const.Planck_constant / cmiext.const.speed_of_light / 100


def inch2m(val):
    """inch -> m"""
    return val * cmiext.const.inch


def invcm2J(val):
    """cm^{-1} -> Joule"""
    return val * 100 * cmiext.codata.codata["Planck constant"][0] * cmiext.codata.codata["speed of light in vacuum"][0]


def m2eV(val):
    """wavelenght (m) -> eV"""
    return 1/ numpy.array(val) / cmiext.codata.codata["electron volt-inverse meter relationship"][0]


def kV_cm2V_m(val):
    """kV/cm -> V/m"""
    return numpy.array(val) / 1e-5

def V_m2kV_cm(val):
    """V/m -> kV/cm"""
    return numpy.array(val) * 1e-5

def A32CM2_V(val):
    """Å^3 -> C M^2 / V"""
    return val*cmiext.const.vacuum_permittivity*((cmiext.const.Angstrom)**3)*4*cmiext.const.pi


# useful strong-field physics conversions
def sfi_au2eV(p, m=cmiext.const.electron_mass):
    """convert electron momentum in atomic units to kinetic energy

    p = momentum in atomic units
    m = mass in SI units
    """
    p = cmiext.convert.numpy.array(p) * cmiext.const.Planck_constant / (2 * numpy.pi * cmiext.const.Bohr_radius)
    return J2eV(p**2 / (2 * m))



def sfi_Up(wl, I=None, E=None):
    """calculate pondermotive energy

    Need so specify exactly one of I or E.

    wl = wavelength in m
    I = intensity in W/m^2
    E = electric field amplitude in V/m

    See https://en.wikipedia.org/wiki/Ponderomotive_energy dor details.
    """
    if ((None == I) and (None == E)) or ((None != I) and (None != E)):
        raise ValueError('sfi_Up: Must specify exactly one of I or E')
    if None == E:
        E = I2E(I)
    w = 2 * numpy.pi * cmiext.const.speed_of_light / wl
    Up = (cmiext.const.electron_charge * E)**2 / (4 * cmiext.const.electron_mass * w**2)
    return Up
