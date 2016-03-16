# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008,2016 Jochen Küpper <jochen.kuepper@cfel.de>
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


__author__ = "Jochen Küpper <jochen.kuepper@cfel.de>"

import numpy as np


class State:
    """State label of molecule

    Currently only asymmetric top notation, but linear tops naturally fit and symmetric tops are placed in as described
    below.

    Public data:
    - max  Upper bound of any individual quantum number - actually any |qn| must be strictly smaller than max

    For symmetric tops only one K quantum number exists, which is placed in (the natural choice of) Ka or Kc depending
    on the case.
    """

    def __init__(self, J=0, Ka=0, Kc=0, M=0, isomer=0):
        self.max = 1000
        self.__initialize(J, Ka, Kc, M, isomer)
        self.__symtop_sign = 1

    def __initialize(self, J=0, Ka=0, Kc=0, M=0, isomer=0):
        """Store all info and reata a unique ID for the state.

        For symmetric tops, :math:`K` is stored in the natural choice of :math:`K_a` or :math:`K_c`. The state depends
        on the sign of the product of :math:`K\cdot{}M`, which is stored in the decimal place 15 (0-14 being used to
        encode J, Ka, Kc, M, isomer).

        """
        assert ((0 <= J < self.max) and (abs(Ka) < self.max) and (abs(Kc) < self.max) and (0 <= M < self.max)
                and (0 <= isomer < self.max))
        self.__labels = np.array([J, Ka, Kc, M, isomer], dtype=np.int64)
        self.__id = np.uint64(0)
        for i in range(self.__labels.size):
            self.__id += np.uint64(abs(self.__labels[i]) * self.max**i)
        # handle negative sign of symmetric-top K*M
        if Ka * Kc * M < 0:
            self.__symtop_sign = -1
            self.__id += np.uint64(1e15)

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

    def nssw(self, forbidden):
        """Give back nuclear spin weight 0 for nuclear-spin-statistically forbidden rve-states, 1 otherwise"""
        if "Ka" == forbidden and self.Ka() % 2 == 1: return 0
        if "Kb" == forbidden and (self.Ka() + self.Kc()) % 2 == 1: return 0
        if "Kc" == forbidden and self.Kc() % 2 == 1: return 0
        return 1

    def symtop_sign():
        return self.__symtop_sign

    def fromid(self, id):
        """Set quantum-numbers form id"""
        id = np.uint64(id)
        self.__id = id
        self.__labels = np.zeros((5,), dtype=np.int64)
        for i in range(5):
            self.__labels[i] = id % self.max
            id //= self.max
        # handle negative sign of symmetric-top K*M
        self.__symtop_sign = id % 10
        id //= 10
        return self

    def fromhdfname(self, hdfname):
        """Set quantum-numbers form hdf name.

        See hdfname() below for a description of the format.

        .. todo:: Implement symtop-sign usage
        """
        name = hdfname.replace("n","-")
        qn = np.array(name.split("/"))
        J, Ka, Kc, M, iso = qn.tolist()
        self.__initialize(np.int64(J[1:]), np.int64(Ka[1:]), np.int64(Kc[1:]), np.int64(M[1:]), np.int64(iso[1:]))
        return self

    def id(self):
        return self.__id

    def name(self):
        return "%d %d %d %d %d" % self.totuple()

    def hdfname(self):
        """Create HDF5 storage file name of state.

        Prepend '_' to all numbers to make them valid Python identifiers. We split the individual quantum numbers by '/'
        in order to provide subgrouping for faster transversal of the HDF5 directory.

        .. todo:: Implement symtop-sign usage
        """
        name = "_%d/_%d/_%d/_%d/_%d" % self.totuple()
        name.replace("-","n")
        return name.replace("-","n")

    def toarray(self):
        return self.__labels

    def tolist(self):
        return self.__labels.tolist()

    def totuple(self):
        return tuple(self.__labels.tolist())



class RVState:
    """State label of molecule for rovibrational states
        
        Currently only asymmetric top notation, but linear tops naturally fit and symmetric tops are placed in as described
        below.
        
        Public data:
        - max  Upper bound of any individual quantum number - actually any |qn| must be strictly smaller than max
        
        For symmetric tops only one K quantum number exists, which is placed in (the natural choice of) Ka or Kc depending
        on the case.
        """
    
    def __init__(self, J=0, Ka=0, Kc=0, M=0, V=0, isomer=0):
        self.max = 1000
        self.__initialize(J, Ka, Kc, M, V, isomer)
    
    def __initialize(self, J=0, Ka=0, Kc=0, M=0, V=0, isomer=0):
        """Store all info and reata a unique ID for the state.
            
            the decimal places 0-17 being used to encode J, Ka, Kc, M, V, isomer.
            """
        assert ((0 <= J < self.max) and (abs(Ka) < self.max) and (abs(Kc) < self.max) and (0 <= M < self.max)
                and (0 <= V <= self.max) and (0 <= isomer < self.max))
        self.__labels = num.array([J, Ka, Kc, M, V, isomer], dtype=num.int64)
        self.__id = num.uint64(0)
        for i in range(self.__labels.size):
            self.__id += num.uint64(abs(self.__labels[i]) * self.max**i)
    
    def J(self):
        return self.__labels[0]
    
    def Ka(self):
        return self.__labels[1]
    
    def Kc(self):
        return self.__labels[2]
    
    def M(self):
        return self.__labels[3]

    def V(self):
        return self.__labels[4]    
    
    def isomer(self):
        return self.__labels[5]
    
    def nssw(self, forbidden):
        """Give back nuclear spin weight 0 for nuclear-spin-statistically forbidden rve-states, 1 otherwise"""
        if "Ka" == forbidden and self.Ka() % 2 == 1: return 0
        if "Kb" == forbidden and (self.Ka() + self.Kc()) % 2 == 1: return 0
        if "Kc" == forbidden and self.Kc() % 2 == 1: return 0
        return 1
    
    def fromid(self, id):
        """Set quantum-numbers form id"""
        self.__id = num.uint64(id)
        self.__labels = num.zeros((6,), dtype=num.int64)
        for i in range(6):
            self.__labels[i] = id % self.max
            id //= self.max
        return self
    
    def fromhdfname(self, hdfname):
        """Set quantum-numbers form hdf name.
            
            See hdfname() below for a description of the format.
            """
        name = hdfname.replace("n","-")
        qn = num.array(name.split("/"))
        J, Ka, Kc, M, V, iso = qn.tolist()
        self.__initialize(num.int64(J[1:]), num.int64(Ka[1:]), num.int64(Kc[1:]), num.int64(M[1:]), num.int64(V[1:]), num.int64(iso[1:]))
        return self
    
    def id(self):
        return self.__id
    
    def name(self):
        return "%d %d %d %d %d %d" % self.totuple()
    
    def hdfname(self):
        """Create HDF5 storage file name of state.
            
            Prepend '_' to all numbers to make them valid Python identifiers. We split the individual quantum numbers by '/'
            in order to provide subgrouping for faster transversal of the HDF5 directory.
            """
        name = "_%d/_%d/_%d/_%d/_%d/_%d" % self.totuple()
        name.replace("-","n")
        return name.replace("-","n")
    
    def toarray(self):
        return self.__labels
    
    def tolist(self):
        return self.__labels.tolist()
    
    def totuple(self):
        return tuple(self.__labels.tolist())



def create_statelist(Jmax):
    """Create a list of |State|s of all asymmetric rotor quantum states up to Jmax.

    The list is not ordered energetically.
    """
    return


# state = State(1, 0, 1, 1)
