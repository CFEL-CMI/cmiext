#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2009 Jochen Küpper <software@jochen-kuepper.de>
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this programm for scientific work, you must correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.
from __future__ import division

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"
__doc__ = """Implementation of Wigner 3j, 6j, 9j symbols and reduced rotation matrixes.

Based on the wigner GSL extension module by Jonathan Underwood.
"""


import ctypes
import ctypes.util

cwigner_filename = ctypes.util.find_library("jkext._wigner.so")
cwigner = ctypes.CDLL(cwigner_filename)


def wigner_drot (j, m1, m2, theta):
    """Reduced Wigner rotation matrix."""
    result = ctypes.c_double()
    status = cwigner.gsl_sf_wigner_drot_e(j, m1, m2, ctypes.c_double(theta), ctypes.byref(result))
    assert 0 == status, """Error in gsl_sf_wigner_drot_e C function."""
    return result
