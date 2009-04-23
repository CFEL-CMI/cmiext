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


class State:
    """State label of molecule (currently only asymmetric top notation)

    Public data:
    - max  Upper bound of any individual quantum number - actually any qn must be strictly smaller than max
    """

    def __init__(self, J=0, Ka=0, Kc=0, M=0, isomer=0):
        self.max = 1000
        self.__initialize(J, Ka, Kc, M, isomer)

    def __initialize(self, J=0, Ka=0, Kc=0, M=0, isomer=0):
        assert (J < self.max) and (Ka < self.max) and (Kc < self.max) and (M < self.max) and (isomer < self.max)
        self.__labels = num.array([J, Ka, Kc, M, isomer], dtype=num.uint64)
        self.__id = num.uint64(0)
        for i in range(self.__labels.size):
            self.__id += num.uint64(self.__labels[i] * self.max**i)

    def J(self):
        return self.__labels[0]

    def Ka(self):
        return self.__labels[1]

    def Kc(self):
        return self.__labels[2]

    def M(self):
        return self.__labels[3]

    def isomer(self):
        return self.__labels[4]

    def fromid(self, id):
        """Set quantum-numbers form id"""
        self.__id = num.uint64(id)
        self.__labels = num.zeros((5,), dtype=num.uint64)
        for i in range(5):
            self.__labels[i] = id % self.max
            id //= self.max
        return self

    def fromhdfname(self, hdfname):
        """Set quantum-numbers form hdf name.

        See hdfname() below for a description of the format.
        """
        qn = num.array(hdfname.split('/'))
        J, Ka, Kc, M, iso = qn.tolist()
        self.__initialize(num.uint64(J[1:]), num.uint64(Ka[1:]), num.uint64(Kc[1:]), num.uint64(M[1:]), num.uint64(iso[1:]))
        return self

    def id(self):
        return self.__id

    def name(self):
        return "%d %d %d %d %d" % self.totuple()

    def hdfname(self):
        """Create HDF5 storage file name of state.

        Prepend '_' to all numbers to make them valid Python identifiers. We split the individual quantum numbers by '/'
        in order to provide subgrouping for faster transversal of the HDF5 directory.
        """
        return "_%d/_%d/_%d/_%d/_%d" % self.totuple()

    def toarray(self):
        return self.__labels

    def tolist(self):
        return self.__labels.tolist()

    def totuple(self):
        return tuple(self.__labels.tolist())
