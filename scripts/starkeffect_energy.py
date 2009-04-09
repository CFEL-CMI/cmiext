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
#
# some SGE commands for batch usage convenience
#$ -cwd
#$ -e $JOB_NAME.sge
#$ -o $JOB_NAME.sge
#$ -S $HOME/.python/bin/python
#$ -V
from __future__ import division

"""Demonstration script to print Stark curves from HDF file to commandline

Copyright (C) 2009 Jochen Küpper"""

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num
import sys

import jkext.molecule as molecule
from jkext.convert import *
from jkext.state import State

# create Molecule object and specify storage file
storagename = "benzonitrile.hdf"
mol = molecule.Molecule(storage=storagename)

for J in range(0, 4):
    Ka = 0
    for Kc in range(J, -1, -1):
        state = State(J, Ka, Kc, 0, 0)
        fields, energies = mol.starkeffect(state)
        print state.name(), V_m2kV_cm(fields), J2Hz(energies) / 1e6
        if Kc > 0:
            Ka += 1
            state = State(J, Ka, Kc, 0, 0)
            fields, energies = mol.starkeffect(state)
            print state.name(), V_m2kV_cm(fields), J2Hz(energies) / 1e6
