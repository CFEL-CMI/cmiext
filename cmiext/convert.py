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

import numpy as np
import scipy
import scipy.constants



def D2Cm(val):
    """Convert dipole moment from Debye to Coulomb * meter"""
    return np.array(val) * 1e-21 / scipy.constants.speed_of_light


def Cm2D(val):
    """Convert dipole moment from Coulomb * meter to Debye"""
    return np.array(val) / (1e-21 / scipy.constants.speed_of_light)


def eV2m(val):
    """eV -> wavelenght (m)"""
    return 1 / np.array(val) / scipy.constants.value('electron volt-inverse meter relationship')


def eV2invcm(val):
    """eV -> wavenumber (cm^{-1})"""
    return 0.01 * np.array(val) * scipy.constants.value('electron volt-inverse meter relationship')


def E2I(val):
    """field amplitude (V/m) -> intensity (m)"""
    return 0.5 * scipy.constants.speed_of_light * scipy.constants.epsilon_0 * np.array(val)**2


def Hz2J(val):
    """Hertz -> Joule"""
    return np.array(val) * scipy.constants.Planck


def MHz2J(val):
    """Mega-Hertz -> Joule"""
    return Hz2J(val * 1e6)


def I2E(val):
    """intensity (W/m^2) -> amplitude (V/m)"""
    return np.sqrt(2 * val / (scipy.constants.speed_of_light * scipy.constants.epsilon_0))


def invcm2Hz(val):
    """wavenumber (cm^{-1}) -> frequency (Hz)"""
    return val * 100 * scipy.constants.speed_of_light


def invcm2J(val):
    """cm^{-1} -> Joule"""
    return invcm2Hz(val) * scipy.constants.Planck


def invcm2m(val):
    """wavenumber (cm^{-1}) -> wavelength (m)"""
    return 0.01 / val


def J2eV(val):
    """Joule -> electron volt"""
    return np.array(val) / scipy.constants.value('electron volt-joule relationship')


def J2Hz(val):
    """Joule -> Hertz"""
    return np.array(val) / scipy.constants.Planck


def J2MHz(val):
    """Joule -> Mega-Hertz"""
    return J2Hz(val) / 1e6


def J2invcm(val):
    """Joule -> cm^{-1}"""
    return val / scipy.constants.Planck / scipy.constants.speed_of_light / 100


def inch2m(val):
    """inch -> m"""
    return val * scipy.constants.inch


def invcm2J(val):
    """cm^{-1} -> Joule"""
    return val * 100 * scipy.constants.Planck * scipy.constants.speed_of_light


def m2eV(val):
    """wavelenght (m) -> eV"""
    return 1/ np.array(val) / scipy.constants.value('electron volt-inverse meter relationship')


def kV_cm2V_m(val):
    """kV/cm -> V/m"""
    return np.array(val) / 1e-5

def V_m2kV_cm(val):
    """V/m -> kV/cm"""
    return np.array(val) * 1e-5

def A32CM2_V(val):
    """polarizability volume (Å^3) -> polarizability (C M^2 / V)"""
    return val * scipy.constants.epsilon_0 * (scipy.constants.angstrom)**3 * 4*scipy.constants.pi


# useful strong-field physics conversions
def sfi_au2eV(p, m=scipy.constants.electron_mass):
    """convert electron momentum in atomic units to kinetic energy in eV

    p = momentum in atomic units
    m = mass in SI units
    """
    p = np.array(p) * scipy.constants.Planck / (2 * np.pi * scipy.constants.value('Bohr radius'))
    return J2eV(p**2 / (2 * m))


# useful strong-field physics conversions
def sfi_Keldysh(Ip, Up):
    """calulate Keldysh parameter from Ip and Up

    Ip = ionization potential
    Up = pondermotive energy
    """
    return scipy.sqrt(Ip/(2*Up))



def sfi_Up(wl, I=None, E=None):
    """calculate pondermotive energy in J(!)

    Need so specify exactly one of I or E.

    wl = wavelength in m
    I = intensity in W/m^2
    E = electric field amplitude in V/m

    See https://en.wikipedia.org/wiki/Ponderomotive_energy for details.
    """
    if ((None == I) and (None == E)) or ((None != I) and (None != E)):
        raise ValueError('sfi_Up: Must specify exactly one of I or E')
    if None == E:
        E = I2E(I)
    w = 2 * np.pi * scipy.constants.speed_of_light / wl
    Up = (scipy.constants.elementary_charge * E)**2 / (4 * scipy.constants.electron_mass * w**2)
    return Up



def sfi_Up_eV(wl, I=None, E=None):
    return cmiext.convert.J2eV(sfi_Up(wl, I, E))
