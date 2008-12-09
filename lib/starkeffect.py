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
from __future__ import division

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num
import numpy.linalg
import convert
from jkext.state import State


class CalculationParameter:
    """Container of parameters for calculation of Stark energies.

    Calculate energy for the specified |fields| (V/m) and rotor type; all calculations are performed in representation
    Ir (x, y, z -> b, c, a).

    General parameters:
    - type: specify the type of rotor, currently only 'A' is implemented.
      - 'A': asymmetric top

    The following parameter are used for an asymmetric top:
    - M, Jmin, Jmax_calc, Jmax_save
    - isomer
    - rotcon (Joule), quartic (Joule), dipole (Coulomb meter)
    - watson, symmetry

    |watson| specifies which reduction of the centrifugal distortion constants of an asymmetric top shall be used.
    - 'A' for Watson's A reduction
    - 'S' for Watson's S reduction

    |symmetry| defines the remaining symmetry of Hamiltonian for the molecule in a DC field. This is used to disentangle
    the block-diagonalization from a Wang transformation of the Hamiltonina matrix.

    It can be None, 'a', 'b', 'c', or 'V' for full Fourgroup symmetry. The latter can only be correct for zero-field
    calculations.
    """
    type = 'A'
    fields = convert.kV_cm2V_m(num.linspace(0., 100., 2)) # V/m
    M = range(0, 2)
    Jmin = 0
    Jmax_calc = 5
    Jmax_save = 2
    isomer = 0
    # molecular parameters
    rotcon = num.zeros((3,), num.float64)    # Joule
    quartic = num.zeros((5,), num.float64)   # Joule
    dipole = num.zeros((3,), num.float64)    # Coulomb meter
    watson=None
    symmetry=None



class AsymmetricRotor:
    """Representation of an asymmetric top for energy level calculation purposes.

    This object will caclulate rotational energies at the specified field strength for the given M-value and J-range and
    all K's.
    """

    def __init__(self, param, M, field):
        """Save the relevant parameters"""
        assert 'A' == param.type.upper()
        # we have not yet calculated the correct energies - mark invalid
        self.__valid = False
        # save parameters internally
        self.__field = num.float64(field)
        self.__rotcon = num.array(param.rotcon, num.float64)
        self.__quartic = num.array(param.quartic, num.float64)
        self.__dipole = num.array(param.dipole, num.float64)
        self.__watson = param.watson
        self.__symmetry = param.symmetry
        # save quantum numbers
        self.__M = int(M) # use the single spefied M
        self.__isomer = int(param.isomer)
        self.__Jmin = int(max(self.__M, param.Jmin))
        self.__Jmax = int(param.Jmax_calc)
        self.__Jmax_save = int(param.Jmax_save)
        self.__matrixsize_Jmin = self.__Jmin *(self.__Jmin-1) + self.__Jmin
        self.__matrixsize = (self.__Jmax + 1) * self.__Jmax + self.__Jmax + 1 - self.__matrixsize_Jmin
        # more checks
        assert self.__rotcon.shape == (3,)
        assert self.__quartic.shape == (5,)
        # some useful constants
        self.__tiny = num.finfo(num.dtype(num.float64)).tiny * 10
        self.__dipole_components = [self.__tiny < abs(self.__dipole[0]),
                                    self.__tiny < abs(self.__dipole[1]),
                                    self.__tiny < abs(self.__dipole[2])]



    def energy(self, state):
        """Return Stark energy for |state|."""
        if self.__valid == False:
            self.__recalculate()
        return self.__levels[state.id()]


    def field(self):
        """Return field for which the Stark energies were calculated."""
        return self.__field


    def states(self):
        """Return list of states for which the Stark energies were calculated."""
        list = []
        M = self.__M
        iso = self.__isomer
        for J in range(self.__Jmin, self.__Jmax_save+1):
            Ka = 0
            for Kc in range(J, -1, -1):
                list.append(State(J, Ka, Kc, M, iso))
                if Kc > 0:
                    Ka += 1
                    list.append(State(J, Ka, Kc, M, iso))
        return list


    def __index(self, J, K):
        blockstart = J*(J-1) + J - self.__matrixsize_Jmin
        return blockstart + K + J


    def __recalculate(self):
        """Needs to be finished!"""
        self.__levels = {}
        self.__full_hamiltonian()
        blocks = self.__wang()
        for symmetry in blocks.keys():
            eval = num.linalg.eigvalsh(blocks[symmetry]) # calculate only energies
            eval = num.sort(eval)
            i = 0
            for state in self.__stateorder(symmetry):
                if state.J() <= self.__Jmax_save:
                    self.__levels[state.id()] = eval[i]
                i += 1
        # done - data is now valid
        self.__valid == True


    def __full_hamiltonian(self):
        # create hamiltonian matrix
        if True == self.__dipole_components[1]: # µ_b != 0
            self.__hmat = num.zeros((self.__matrixsize, self.__matrixsize), num.complex128)
        else:
            self.__hmat = num.zeros((self.__matrixsize, self.__matrixsize), num.float64)
        # start matrix with appropriate field-free rigid-rotor terms
        self.__rigid()
        # add appropriate field-free centrifugal distortion terms
        if self.__watson == 'A':
            self.__watson_A()
        elif self.__watson == 'S':
            self.__watson_S()
        else:
            assert self.__watson == None
        # fill matrix with appropriate Stark terms for nonzero fields
        if self.__tiny < abs(self.__field):
            self.__stark()


    def __rigid(self):
        """Add the rigid-rotor matrix element terms to self.__hmat"""
        sqrt = num.sqrt
        A, B, C = self.__rotcon.tolist()
        for J in range(self.__Jmin, self.__Jmax+1):
            for K in range(-J, J+1):
                self.__hmat[self.__index(J, K), self.__index(J, K)] += (B+C)/2 * (J*(J+1) - K**2) + A * K**2
            for K in range (-J, J-2+1):
                value = (B-C)/4 * sqrt(J*(J+1) - (K+1)*(K+2)) * sqrt(J*(J+1) - K*(K+1))
                self.__hmat[self.__index(J, K+2), self.__index(J, K)] += value
                self.__hmat[self.__index(J, K), self.__index(J, K+2)] += value


    def __stark(self):
        """Add the Stark-effect matrix element terms to self.__hmat"""
        sqrt = num.sqrt
        field = self.__field
        M = self.__M
        muA, muB, muC = self.__dipole
        if self.__dipole_components[0]:
            # matrix elements involving µ_a
            for J in range(self.__Jmin, self.__Jmax):
                for K in range(-J, J+1):
                    if 0 != J:
                        self.__hmat[self.__index(J, K), self.__index(J, K)] += -muA * field * M * K / (J*(J+1))
                        value = (-muA * field * sqrt((J+1)**2 - K**2) * sqrt((J+1)**2 - M**2)
                                  / ((J+1) * sqrt((2*J+1) * (2*J+3))))
                        self.__hmat[self.__index(J+1, K), self.__index(J, K)] += value
                        self.__hmat[self.__index(J, K), self.__index(J+1, K)] += value
            # final diagonal elements
            J = self.__Jmax
            for K in range(-J, J+1):
                self.__hmat[self.__index(J, K), self.__index(J, K)] += -1. * M * K / (J*(J+1)) * muA * field
        if self.__dipole_components[2]:
            # matrix elements involving µ_c
            for J in range(self.__Jmin, self.__Jmax):
                for K in range(-J, J+1):
                    if 0 != J:
                        value = -1. * M * muC * field * (sqrt((J-K) * (J+K+1) ) ) / ( 2*J*(J+1))
                        self.__hmat[self.__index(J, K+1), self.__index(J, K)] += value
                        self.__hmat[self.__index(J, K), self.__index(J, K+1)] += value
                        #K+1 case
                        value = (sqrt( (J+K+1) * (J+K+2) ) * sqrt( (J+1)*(J+1) - M**2  )
                                 / ( 2*(J+1) * sqrt( (2*J+1) * (2*J+3) ) ) * muC * field)
                        self.__hmat[self.__index(J+1, K+1), self.__index(J, K)] += value
                        self.__hmat[self.__index(J, K), self.__index(J+1, K+1)] += value
                        #K-1 case
                        value = ((-1)*sqrt( (J-K+1) * (J-K+2) ) * sqrt( (J+1)**2 - M**2  )
                                 / ( 2*(J+1) * sqrt( (2*J+1) * (2*J+3) ) ) * muC * field)
                        self.__hmat[self.__index(J+1, K-1), self.__index(J, K)] += value
                        self.__hmat[self.__index(J, K), self.__index(J+1, K-1)] += value
        if  self.__dipole_components[1]:
            # matrix elements involving µ_b
            for J in range(self.__Jmin, self.__Jmax):
                for K in range(-J, J+1):
                    if 0 != J:
                        value = 1j* M * muB * field * sqrt((J-K) * (J+K+1)) / (2*J*(J+1))
                        self.__hmat[self.__index(J, K+1), self.__index(J, K)] += value
                        self.__hmat[self.__index(J, K), self.__index(J, K+1)] += value
                        # J+1, K+1 / J-1, K-1 case
                        value = (-1j * muB * field * sqrt((J+K+1) * (J+K+2)) * sqrt((J+1)**2 - M**2)
                                  / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                        self.__hmat[self.__index(J+1, K+1), self.__index(J, K)] += value
                        self.__hmat[self.__index(J, K), self.__index(J+1, K+1)] += value
                        # J+1, K-1 / J-1, K+1 case
                        value = (-1j  * muB * field * sqrt((J-K+1) * (J-K+2)) * sqrt((J+1)**2 - M**2)
                                  / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                        self.__hmat[self.__index(J+1, K-1), self.__index(J, K)] += value
                        self.__hmat[self.__index(J, K), self.__index(J+1, K-1)] += value


    def __stateorder(self, symmetry):
        """Return a list with all states for the given |symmetry| and the current calculation parameters (Jmin, Jmax).

        Needs to be finished!
        """
        return self.states()


    def __wang(self):
        """Wang transform matrix and return a dictionary with the individual (sub)matrices."""
        blocks = {}
        # set up Wang matrix
        Wmat = num.zeros(self.__hmat.shape, num.float64)
        value = 1/num.sqrt(2.)
        for J in range(self.__Jmin, self.__Jmax + 1):
            for K in range(-J, 0):
                Wmat[self.__index(J,  K), self.__index(J,  K)] = -value
                Wmat[self.__index(J, -K), self.__index(J,  K)] = value
                Wmat[self.__index(J,  K), self.__index(J, -K)] = value
                Wmat[self.__index(J, -K), self.__index(J, -K)] = value
            Wmat[self.__index(J, 0), self.__index(J, 0)] = 1.
        # transform Hamiltonian matrix
        self.__hmat = num.dot(num.dot(Wmat, self.__hmat), Wmat)
        # sort out matrix blocks
        if self.__symmetry == None:
            # nothing to do, return
            blocks['N'] = self.__hmat
        elif self.__symmetry == 'V':
            # full Fourgroup symmetry
            #
            # I^r representation, Wang transformed Hamiltonian factorizes into four submatrices E-, E+, O-, O+
            blocks['Ep'] = self.__hmat
            blocks['Em'] = self.__hmat
            blocks['Op'] = self.__hmat
            blocks['Om'] = self.__hmat
        elif self.__symmetry == 'a':
            # C2 rotation about a-axis is symmetry element
            #
            # I^r representation, Wang transformed Hamiltonian factorizes into two submatrices E (contains E- and
            # E+) and O (contains O- and O+). In this case E and O corresponds to columns with Ka even and odd,
            # respectively.
            blocks['E'] = self.__hmat
            blocks['O'] = self.__hmat
        elif self.__symmetry == 'b': # C2 rotation about b-axis is symmetry element
            raise NotImplementedError("Hamiltonian symmetry 'b' not implemented yet")
        elif self.__symmetry == 'c': # C2 rotation about c-axis is symmetry element
            raise NotImplementedError("Hamiltonian symmetry 'c' not implemented yet")
        else:
            # something went wrong
            raise "unknown Hamiltonian symmetry"
        return blocks


    def __watson_A(self):
        """Add the centrifugal distortion matrix element terms in Watson's A reduction to self.__hmat."""
        sqrt = num.sqrt
        DJ, DJK, DK, dJ, dK = self.__quartic.tolist()
        for J in range(self.__Jmin, self.__Jmax+1):
            for K in range(-J, J+1):
                value = DJ * (J*(J+1))**2 - DJK * J*(J+1) * K**2 - DK * K**4
                self.__hmat[self.__index(J, K), self.__index(J, K)] += value
            for K in range (-J, J-2+1):
                value = ((-dJ * J*(J+1) - dK/2 * (K+2)**2 * K**2)
                         * sqrt(J*(J+1) - (K+1)*(K+2)) * sqrt(J*(J+1) - K*(K+1)))
                self.__hmat[self.__index(J, K+2), self.__index(J, K)] += value
                self.__hmat[self.__index(J, K), self.__index(J, K+2)] += value


    def __watson_S(self):
        """Add the centrifugal distortion matrix element terms in Watson's S reduction to self.__hmat."""
        raise NotImplementedError("Watson's S-reduction is not implemented (yet)")



# some simple tests
if __name__ == "__main__":
    print
    p = CalculationParameter
    p.rotcon = num.array([5e9, 2e9, 1.8e9])
    top = AsymmetricRotor(p, 0, 0)
    print top.energy(State(1, 0, 1, 0, 0)) / 1e6
