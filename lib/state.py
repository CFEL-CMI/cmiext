#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2008 Jochen Küpper <software@jochen-kuepper.de>
# see LICENSE file for details

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num


class State:
    """State label of molecule (currently only asymmetric top notation)

    Public data:
    - max  Upper bound of any individual quantum number
    """

    def __init__(self, J=0, Ka=0, Kc=0, M=0, isomer=0):
        self.__labels = num.array([J, Ka, Kc, M, isomer], dtype=num.uint64)
        self.max = 1000
        assert J < self.max and Ka < self.max and Kc < self.max and Kc < self.max and isomer < self.max
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

    def id(self):
        return self.__id

    def name(self):
        return "%d %d %d %d %d" % self.totuple()

    def hdfname(self):
        return "_%d_%d_%d_%d_%d_" % self.totuple()

    def toarray(self):
        return self.__labels

    def tolist(self):
        return self.__labels.tolist()

    def totuple(self):
        return tuple(self.__labels.tolist())



### Local Variables:
### fill-column: 132
### End:
