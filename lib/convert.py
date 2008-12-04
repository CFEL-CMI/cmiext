# -*- coding: utf-8 -*-
#
# Copyright (C) 2008 Jochen Küpper <software@jochen-kuepper.de>
# see LICENSE file for details

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



### Local Variables:
### fill-column: 132
### End:
