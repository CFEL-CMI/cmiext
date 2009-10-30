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
import string


def readVLArray(file, name):
    array = file.getNode(name)
    return num.array(array.read())[0]

def readVLArraysfromsubgroups(file, name, array_name):
    """
    read all arrays of the same type from subgroups of the specified group also return the dirs converted to
    floats assuming they are acfields
    """ 
    groups = file.getNode(name)._v_groups
    acfields = [];
    for acfieldstr in groups:
        acfield = float(string.replace(acfieldstr,'d','.')[1:])
        acfields.append(acfield)
        for array in file.listNodes(name + "/" + acfieldstr,classname='VLArray'):
            if array.name == array_name:
                try:
                    stackedarrays = num.vstack([stackedarrays,array])
                except UnboundLocalError:
                    stackedarrays = array
    stackedarrays[:] = stackedarrays[num.argsort(acfields)]
    acfields = num.sort(acfields)
    return stackedarrays,acfields


def readSubdirs(file, name):
    """
    read all arrays of the same type from subgroups of the specified group also return the dirs converted to
    floats assuming they are acfields
    """ 
    groups = file.getNode(name)._v_groups
    acfields = [];
    for acfieldstr in groups:
        acfield = float(string.replace(acfieldstr,'d','.')[1:])
        acfields.append(acfield)
    acfields = num.sort(acfields)
    return acfields

def writeVLArray(file, groupname, leafname, data, comment="", atom=tables.Float64Atom(shape=()),
                 filters=tables.Filters(complevel=1, complib='zlib')):
    """
    Write a single array, corresponding to a single Stark curve, to the storage file.

    We only use zlib-compression at level 1, because that's apparently as good as any higher level, but should be
    faster. Moreover, we rely on PyTables automatically turning on the HDF5 shuffle filter, what it does when any
    compression is turned on.
    """
    # make sure the group-tree exists
    group = file.root
    for name in groupname.split('/'):
        try:
            group = file.getNode(group, name)
        except tables.exceptions.NodeError:
            try:
                group = file.createGroup(group, name)
            except tables.exceptions.NodeError:
                assert False, "Stark storage error: cannot create non-existing group %s from %s!" % (name, groupname)
    # if the dataset exists already, delete it
    try:
        file.removeNode(group, leafname)
    except tables.exceptions.NodeError, tables.exceptions.NoSuchNodeError:
        pass
    array = file.createVLArray(group, leafname, atom, comment, filters)
    array.append(data)
    array.flush()
