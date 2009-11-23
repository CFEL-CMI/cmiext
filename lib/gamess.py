# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2009 Jochen Küpper <software@jochen-kuepper.de>
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

import numpy as num
import molecule

def atom(line):
    name, charge, x, y, z = line.split()
    return molecule.Atom(name, (x, y, z))


def geometries(file):
    geom = []
    new_structure = False
    structure = []
    for line in open(file):
        if 0 < line.find("COORDINATES OF ALL ATOMS ARE"):
            new_structure = True
            structure = []
        elif len(line) < 2:
            new_structure = False
            geom.append(structure)
        elif True == new_structure:
            if (0 < line.find("ATOM")) or (0 < line.find("---")):
                continue
            structure.append(atom(line))
    return geom
