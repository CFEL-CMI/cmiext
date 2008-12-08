#!/usr/bin/env python
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
#
# some SGE commands for batch usage convenience
#$ -cwd
#$ -e $JOB_NAME.sge
#$ -o $JOB_NAME.sge
#$ -S $HOME/.python/bin/python
#$ -V
from __future__ import division

"""Stark effect calculations for benzonitrile, based on constants from Wohlfart, Schnell, Grabow,
Küpper, J. Mol. Spec. 247, 119-121 (2008)

Copyright (C) 2008 Jochen Küpper"""

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num
import jkext as jk
import jkext.convert as convert
import jkext.molecule as molecule
import jkext.starkeffect as starkeffect

# create Molecule object and specify storage file
bn = molecule.Molecule(storage="benzonitrile.hdf")

# set molecular parameters
param = starkeffect.CalculationParameter
param.isomer = 0
param.watson = 'A'
param.symmetry = None # 'a'
param.rotcon = convert.Hz2J(num.array([5655.2654e6, 1546.875864e6, 1214.40399]))
param.quartic = convert.Hz2J(num.array([45.6, 938.1, 500, 10.95, 628]))
param.dipole = convert.D2Cm(num.array([4.5152, 0., 0.]))
# calculation details
param.M = [0]
param.Jmin = 0
param.Jmax_calc = 5
param.Jmax_save = 3
param.fields = convert.kV_cm2V_m(num.linspace(0., 100., 2))

# perform calculation
bn.starkeffect_calculation(param)
