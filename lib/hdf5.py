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

import numpy as num
import tables


def readVLArray(file, name):
    array = file.getNode(name)
    return num.array(array.read())[0]


def writeVLArray(file, groupname, leafname, data, comment="", atom=tables.Float64Atom(shape=()),
                 filters=tables.Filters(complevel=9, complib='zlib')):
    root = file.root
    # make sure the group exists
    try:
        group = file.getNode(root, groupname)
    except tables.exceptions.NodeError:
        group = root
        for name in groupname.split('/'):
            if 0 < len(name):
                try:
                    group = file.createGroup(group, name)
                except tables.exceptions.NodeError:
                    pass
    # if the dataset exists already, delete it
    try:
        file.removeNode(group, leafname)
    except tables.exceptions.NodeError, tables.exceptions.NoSuchNodeError:
        pass
    array = file.createVLArray(group, leafname, atom, comment, filters)
    array.append(data)
