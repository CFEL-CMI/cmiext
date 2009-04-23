#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2008,2009 Jochen Küpper <software@jochen-kuepper.de>
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

# really use scipy as numpy, so we are sure we use Fortran codes of eigvalsh and dgemm
import scipy as num
import scipy.linalg.fblas

import jkext.convert
from jkext.state import State


class CalculationParameter:
    """Container of parameters for calculation of Stark energies.

    Calculate energy for the specified |dcfields| (V/m) and rotor type; all calculations are performed in representation
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
    M = range(0, 2)
    Jmax_calc = 5
    Jmax_save = 2
    isomer = 0
    # fields
    acfields = num.zeros((1,), num.float64)
    dcfields = jkext.convert.kV_cm2V_m(num.array((0, 100.), num.float64))
    # molecular parameters
    rotcon = num.zeros((3,), num.float64)    # Joule
    quartic = num.zeros((5,), num.float64)   # Joule
    dipole = num.zeros((3,), num.float64)    # Coulomb meter
    polarizability = num.zeros((3,3), num.float64)
    watson=None
    symmetry='N'



class AsymmetricRotor:
    """Representation of an asymmetric top for energy level calculation purposes.

    This object will caclulate rotational energies at the specified AC and DC field strength for the given M-value and
    J-range and all K's.
    """

    def __init__(self, param, M, acfield=0., dcfield=0.):
        """Save the relevant parameters"""
        assert 'A' == param.type.upper()
        # we have not yet calculated the correct energies - mark invalid
        self.__valid = False
        self.__stateorder_valid = False
        # save parameters internally
        self.__acfield = num.float64(acfield)
        self.__dcfield = num.float64(dcfield)
        self.__rotcon = num.array(param.rotcon, num.float64)
        self.__quartic = num.array(param.quartic, num.float64)
        self.__dipole = num.array(param.dipole, num.float64)
        self.__watson = param.watson
        self.__symmetry = param.symmetry # symmetry of Hamiltonian (possible values: 'N', 'C2a', 'C2b', 'C2c', 'V')
        # save quantum numbers
        self.__M = int(M) # use the single spefied M
        self.__isomer = int(param.isomer)
        self.__Jmin = self.__M # this must be equal to self.__M (in Stark calculation all J couple)
        self.__Jmax = int(param.Jmax_calc)
        self.__Jmax_save = int(param.Jmax_save)
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


    def field_AC(self):
        """Return AC field for which the Stark energies were calculated."""
        return self.__acfield


    def field_DC(self):
        """Return DC field for which the Stark energies were calculated."""
        return self.__dcfield


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
        # this requires a correct "global" value of self.__Jmin_matrixsize, which is set in __full_hamiltonian.
        # Therefore, we must be called only through __full_hamiltonian
        blockstart = J*(J-1) + J - self.__Jmin_matrixsize
        return blockstart + K + J


    def __recalculate(self):
        """Perform calculation of rotational state energies for current parameters"""
        self.__levels = {}
        blocks = self.__full_hamiltonian(self.__Jmin, self.__Jmax, self.__acfield, self.__dcfield, self.__symmetry)
        for symmetry in blocks.keys():
            eval = num.linalg.eigvalsh(blocks[symmetry]) # calculate only energies
            eval = num.sort(eval)
            i = 0
            for state in self.__stateorder(symmetry):
                if state.J() <= self.__Jmax_save:
                    self.__levels[state.id()] = eval[i]
                i += 1
        # done - data is now valid
        self.__valid = True


    def __full_hamiltonian(self, Jmin, Jmax, acfield, dcfield, symmetry):
        """Return block-diagonalized Hamiltonian matrix (blocks)"""
        self.__Jmin_matrixsize = Jmin *(Jmin-1) + Jmin
        matrixsize = (Jmax + 1) * Jmax + Jmax + 1 - self.__Jmin_matrixsize
        # create hamiltonian matrix
        if True == self.__dipole_components[2]:
            # if µ_c != 0, the matrix is complex (and hermitean)
            hmat = num.zeros((matrixsize, matrixsize), num.complex128)
        else:
            # otherwise the matrix is real (and symmetric)
            hmat = num.zeros((matrixsize, matrixsize), num.float64)
        # start matrix with appropriate field-free rigid-rotor terms
        self.__rigid(hmat, Jmin, Jmax)
        # add appropriate field-free centrifugal distortion terms
        if self.__watson == 'A':
            self.__watson_A(hmat, Jmin, Jmax)
        elif self.__watson == 'S':
            self.__watson_S(hmat, Jmin, Jmax)
        else:
            assert self.__watson == None
        # fill matrix with appropriate Stark terms for nonzero fields
        if None != acfield and self.__tiny < abs(acfield):
            self.__stark_AC(hmat, Jmin, Jmax, acfield)
        if None != dcfield and self.__tiny < abs(dcfield):
            self.__stark_DC(hmat, Jmin, Jmax, dcfield)
        blocks = self.__wang(hmat, symmetry, Jmin, Jmax)
        del hmat
        return blocks


    def __rigid(self, hmat, Jmin, Jmax):
        """Add the rigid-rotor matrix element terms to hmat"""
        sqrt = num.sqrt
        A, B, C = self.__rotcon.tolist()
        for J in range(Jmin, Jmax+1):
            for K in range(-J, J+1):
                hmat[self.__index(J, K), self.__index(J, K)] += (B+C)/2 * (J*(J+1) - K**2) + A * K**2
            for K in range (-J, J-2+1):
                value = (B-C)/4 * sqrt((J*(J+1) - K*(K+1)) * (J*(J+1) - (K+1)*(K+2)))
                hmat[self.__index(J, K+2), self.__index(J, K)] += value
                hmat[self.__index(J, K), self.__index(J, K+2)] += value


    def __stark_AC(self, hmat, Jmin, Jmax, acfield):
        """Add the ac Stark-effect matrix element terms to hmat"""
        pass


    def __stark_DC(self, hmat, Jmin, Jmax, dcfield):
        """Add the dc Stark-effect matrix element terms to hmat"""
        sqrt = num.sqrt
        M = self.__M
        muA, muB, muC = self.__dipole
        if self.__dipole_components[0]:
            # matrix elements involving µ_a
            for J in range(Jmin, Jmax):
                for K in range(-J, J+1):
                    if 0 != J:
                        hmat[self.__index(J, K), self.__index(J, K)] += -muA * dcfield * M * K / (J*(J+1))
                    value = (-muA * dcfield * sqrt((J+1)**2 - K**2) * sqrt((J+1)**2 - M**2)
                              / ((J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.__index(J+1, K), self.__index(J, K)] += value
                    hmat[self.__index(J, K), self.__index(J+1, K)] += value
            # final diagonal elements
            J = Jmax
            for K in range(-J, J+1):
                hmat[self.__index(J, K), self.__index(J, K)] += -1. * M * K / (J*(J+1)) * muA * dcfield
        if self.__dipole_components[1]:
            # matrix elements involving µ_b
            for J in range(Jmin, Jmax):
                for K in range(-J, J+1):
                    if 0 != J:
                        value = -1 * M * muB * dcfield * (sqrt((J-K) * (J+K+1) ) ) / (2*J*(J+1))
                        hmat[self.__index(J, K+1), self.__index(J, K)] += value
                        hmat[self.__index(J, K), self.__index(J, K+1)] += value
                    # J+1, K+1 / J-1, K-1 case
                    value = (muB * dcfield * sqrt(((J+K+1) * (J+K+2)) * ((J+1)**2 - M**2))
                             / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.__index(J+1, K+1), self.__index(J, K)] += value
                    hmat[self.__index(J, K), self.__index(J+1, K+1)] += value
                    # J+1, K-1 / J-1, K+1 case
                    value = (-1 * muB * dcfield * sqrt(((J-K+1) * (J-K+2)) * ((J+1)**2 - M**2))
                              / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.__index(J+1, K-1), self.__index(J, K)] += value
                    hmat[self.__index(J, K), self.__index(J+1, K-1)] += value
        if  self.__dipole_components[2]:
            # matrix elements involving µ_c
            for J in range(Jmin, Jmax):
                for K in range(-J, J+1):
                    if 0 != J:
                        value = 1j* M * muC * dcfield * sqrt((J-K) * (J+K+1)) / (2*J*(J+1))
                        hmat[self.__index(J, K+1), self.__index(J, K)] += value
                        hmat[self.__index(J, K), self.__index(J, K+1)] += value
                    # J+1, K+1 / J-1, K-1 case
                    value = (-1j * muC * dcfield * sqrt((J+K+1) * (J+K+2)) * sqrt((J+1)**2 - M**2)
                              / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.__index(J+1, K+1), self.__index(J, K)] += value
                    hmat[self.__index(J, K), self.__index(J+1, K+1)] += value
                    # J+1, K-1 / J-1, K+1 case
                    value = (-1j  * muC * dcfield * sqrt((J-K+1) * (J-K+2)) * sqrt((J+1)**2 - M**2)
                              / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.__index(J+1, K-1), self.__index(J, K)] += value
                    hmat[self.__index(J, K), self.__index(J+1, K-1)] += value


    def __stateorder(self, symmetry):
        """Return a list with all states for the given |symmetry| and the current calculation parameters (Jmin, Jmax).

        Needs to be finished for C2b and C2c! Split labels and eigenvlaues into eight blocks for each of J Ka Kc 'e' or
        'o'. Then we can merge them differently for J even and odd, such that they provide the correct state labeling.
        """
        if False == self.__stateorder_valid:
            self.__stateorder_dict = {}
            M = self.__M
            iso = self.__isomer
            eigenvalues = {'A': [], 'Ba': [], 'Bb': [], 'Bc': []}
            label = {'A': [], 'Ba': [], 'Bb': [], 'Bc': []}
            for J in range(M, self.__Jmax+1):
                Ka = 0
                for Kc in range(J,-1,-1):
                    if Ka%2 == 0 and Kc%2 == 0:  label['A'].append(State(J, Ka, Kc, M, iso))  # ee
                    elif Ka%2 == 0 and Kc%2 !=0: label['Ba'].append(State(J, Ka, Kc, M, iso)) # eo
                    elif Ka%2 != 0 and Kc%2 !=0: label['Bb'].append(State(J, Ka, Kc, M, iso)) # oo
                    else:                        label['Bc'].append(State(J, Ka, Kc, M, iso)) # oe (Ka%2 != 0 and Kc%2 == 0)
                    if Kc > 0:
                        Ka = Ka+1
                        if Ka%2 == 0 and Kc%2 == 0:   label['A'].append(State(J, Ka, Kc, M, iso))  # ee
                        elif Ka%2 == 0 and Kc%2 != 0: label['Ba'].append(State(J, Ka, Kc, M, iso)) # eo
                        elif Ka%2 != 0 and Kc%2 != 0: label['Bb'].append(State(J, Ka, Kc, M, iso)) # oo
                        else:                         label['Bc'].append(State(J, Ka, Kc, M, iso)) # oe (Ka%2 != 0 and Kc%2 == 0)
                # get block diagonal hamiltonian (make sure you calculate this in 'V'!)
                if 0 == J:
                    blocks = {'A': num.zeros((1, 1), num.float64)}
                else:
                    blocks = self.__full_hamiltonian(J, J, None, None, 'V')
                # store sorted eigenenergies for respective J and block
                for sym in blocks.keys():
                    if 0 < blocks[sym].size:
                        eigenvalues[sym] += num.sort(num.linalg.eigvalsh(num.array(blocks[sym]))).tolist()
            for sym in label.keys():
                print "sym %2s:" % (sym,),
                for state in label[sym]:
                    print state.name(), ";",
                print
            # sort assignments according to energy
            if 'V' == self.__symmetry:
                symmetries = ['A', 'Ba', 'Bb', 'Bc']
            elif 'C2a' == self.__symmetry:
                eigenvalues['Aa'] = eigenvalues['A'] + eigenvalues['Ba']
                eigenvalues['bc'] = eigenvalues['Bb'] + eigenvalues['Bc']
                label['Aa'] = label['A'] + label['Ba']
                label['bc'] = label['Bb'] + label['Bc']
                symmetries = ['Aa', 'bc']
            elif 'C2b' == self.__symmetry:
                eigenvalues['Ab'] = eigenvalues['A'] + eigenvalues['Bb']
                eigenvalues['ac'] = eigenvalues['Bb'] + eigenvalues['Bc']
                label['Ab'] = label['A'] + label['Ba']
                label['ac'] = label['Bb'] + label['Bc']
                symmetries = ['Ab', 'ac']
            elif 'C2c' == self.__symmetry:
                eigenvalues['Ac'] = eigenvalues['A'] + eigenvalues['Bc']
                eigenvalues['ab'] = eigenvalues['Ba'] + eigenvalues['Bb']
                label['Ac'] = label['A'] + label['Bc']
                label['ab'] = label['Bb'] + label['Bb']
                symmetries = ['Ac', 'ab']
            elif 'N' == self.__symmetry:
                eigenvalues['N'] = eigenvalues['A'] + eigenvalues['Ba'] + eigenvalues['Bb'] + eigenvalues['Bc']
                label['N'] = label['A'] + label['Ba'] + label['Bb'] + label['Bc']
                symmetries = ['N']
            else:
                raise NotImplementedError("Hamiltonian symmetry %s not implemented" % (self.__symmetry, ))
            if 'V' != self.__symmetry:
                # free unused memories
                del label['A'], label['Ba'], label['Bb'], label['Bc']
            for sym in symmetries:
                idx = num.argsort(eigenvalues[sym])
                self.__stateorder_dict[sym] = num.array(label[sym])[idx]
            self.__stateorder_valid = True
        return self.__stateorder_dict[symmetry]


    def __wang(self, hmat, symmetry, Jmin, Jmax):
        """Wang transform matrix and return a dictionary with the individual (sub)matrices."""
        matrixsize = ((Jmax + 1) * Jmax + Jmax + 1) - (Jmin *(Jmin-1) + Jmin)
        blocks = {}
        # set up Wang matrix
        Wmat = num.zeros(hmat.shape, num.float64)
        value = 1/num.sqrt(2.)
        for J in range(Jmin, Jmax+1):
            for K in range(-J, 0):
                Wmat[self.__index(J,  K), self.__index(J,  K)] = -value
                Wmat[self.__index(J, -K), self.__index(J,  K)] = value
                Wmat[self.__index(J,  K), self.__index(J, -K)] = value
                Wmat[self.__index(J, -K), self.__index(J, -K)] = value
            Wmat[self.__index(J, 0), self.__index(J, 0)] = 1.
        # transform Hamiltonian matrix
        dot = lambda a, b: scipy.linalg.fblas.dgemm(1., a, b)
        hmat = dot(dot(Wmat, hmat), Wmat)
        print "Wang-transformed Hamiltonian:"
        for i in range(hmat.shape[0]):
            for j in range(hmat.shape[0]):
                print "%10.3g" % (abs(hmat[i,j]), ),
            print("\n")
        # delete Wang matrix (it's not used anymore)
        del Wmat
        # sort out matrix blocks
        if 'N' == symmetry:
            # nothing to do, return
            blocks['N'] = hmat
        elif symmetry == 'V':
            # full Fourgroup symmetry (only field free Hamiltonian)
            # I^r representation, Wang transformed Hamiltonian factorizes into four submatrices E-, E+, O-, O+,
            # or, as used here, A, Ba, Bb, Bc
            idx = {'A': [], 'Ba': [], 'Bb': [], 'Bc': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0
                        if 0 == K % 2: order.append('Ba') # K even
                        else: order.append('Bc') # K odd
                    for K in range(0, J+1): # K >= 0
                        if 0 == K % 2: order.append('A') # K even
                        else: order.append('Bb') # K odd
                else: # J odd
                    for K in range(-J, 0): # K < 0
                        if 0 == K % 2: order.append('A') # K even
                        else: order.append('Bb') # K odd
                    for K in range(0, J+1): # K >= 0
                        if 0 == K % 2: order.append('Ba') # K even
                        else: order.append('Bc') # K odd
                for j in range(2*J+1):
                    idx[order[j]].append(i+j)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        elif symmetry == 'C2a':
            # C2 rotation about a-axis is symmetry element
            #
            # I^r representation, Wang transformed Hamiltonian factorizes into two submatrices E = Aa (contains E+ and
            # E- / A and Ba) and O (contains O+ and O- / Bb And Bc).
            # In this case E and O corresponds to columns with K even and odd, respectively.
            idx = {'Aa': [], 'bc': []}
            if 0 == Jmin % 2: # Jmin even
                order = ['Aa', 'bc']
            else: # J odd
                order = ['bc', 'Aa']
            for i in range(matrixsize):
                idx[order[i%2]].append(i)
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
#             for key in blocks:
#                 print "Wang-transformed block:", key
#                 block = blocks[key]
#                 for i in range(block.shape[0]):
#                     for j in range(block.shape[1]):
#                         print "%10.3g" % (abs(block[i,j]),),
#                     print("\n")
        elif symmetry == 'C2b': # C2 rotation about b-axis is symmetry element
            # C2 rotation about b-axis is symmetry element
            #
            # I^r representation, Wang transformed Hamiltonian factorizes into two submatrices 'Ab' (contains 'A' and 'Bb')
            # and 'ac' (contains 'Ba' and 'Bc').
            idx = {'Ab': [], 'ac': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0): # K < 0
                        if 0 == K % 2: order.append('ac') # K even
                        else: order.append('Ab') # K odd
                    for K in range(0, J+1): # K >= 0
                        if 0 == K % 2: order.append('Ab') # K even
                        else: order.append('ac') # K odd
                else: # J odd
                    for K in range(-J, 0): # K < 0
                        if 0 == K % 2: order.append('Ab') # K even
                        else: order.append('ac') # K odd
                    for K in range(0, J+1): # K >= 0
                        if 0 == K % 2: order.append('ac') # K even
                        else: order.append('Ab') # K odd
                for j in range(2*J+1):
                    idx[order[j]].append(i+j)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
#             for key in blocks:
#                 print "Wang-transformed block:", key
#                 block = blocks[key]
#                 for i in range(block.shape[0]):
#                     for j in range(block.shape[1]):
#                         print "%10.5g" % (block[i,j],),
#                     print("\n")
        elif symmetry == 'C2c':
            # C2 rotation about c-axis is symmetry element
            #
            # I^r representation, Wang transformed Hamiltonian factorizes into two submatrices 'Ac' (contains 'A' and
            # 'Bc') and 'ab' (contains 'Ba' and 'Bb').
            idx = {'Ac': [], 'ab': []}
            i = 0
            for J in range(Jmin, Jmax+1):
                order = []
                if 0 == J % 2: # J even
                    for K in range(-J, 0):  order.append('ab')
                    for K in range(0, J+1): order.append('Ac')
                else: # J odd
                    for K in range(-J, 0):  order.append('Ac')
                    for K in range(0, J+1): order.append('ab')
                for j in range(2*J+1):
                    idx[order[j]].append(i+j)
                i += 2*J+1
            for sym in order:
                if 0 < len(idx[sym]):
                    blocks[sym] = hmat[num.ix_(idx[sym], idx[sym])]
        else:
            # something went wrong
            raise SyntaxError("unknown Hamiltonian symmetry")
        return blocks


    def __watson_A(self, hmat, Jmin, Jmax):
        """Add the centrifugal distortion matrix element terms in Watson's A reduction to hmat."""
        matrixsize_Jmin = Jmin *(Jmin-1) + Jmin
        sqrt = num.sqrt
        DJ, DJK, DK, dJ, dK = self.__quartic.tolist()
        for J in range(Jmin, Jmax+1):
            for K in range(-J, J+1):
                value = -DJ * (J*(J+1))**2 - DJK * J*(J+1)*K**2 - DK * K**4
                hmat[self.__index(J, K), self.__index(J, K)] += value
            for K in range (-J, J-2+1):
                value = ((-dJ * J*(J+1) - dK/2 * ((K+2)**2 + K**2))
                         * sqrt((J*(J+1) - K*(K+1)) * (J*(J+1) - (K+1)*(K+2))))
                hmat[self.__index(J, K+2), self.__index(J, K)] += value
                hmat[self.__index(J, K), self.__index(J, K+2)] += value


    def __watson_S(self):
        """Add the centrifugal distortion matrix element terms in Watson's S reduction to hmat."""
        raise NotImplementedError("Watson's S-reduction is not implemented (yet)")



# some simple tests
if __name__ == "__main__":
    print
    p = CalculationParameter
    p.Jmax_calc =  3
    p.Jmax_save =  3
    p.M = [0]
    p.rotcon = jkext.convert.Hz2J(num.array([5e9, 2e9, 1.4e9]))
    #p.quartic = jkext.convert.Hz2J([1e3, 1e3, 1e3, 1e3, 1e3])
    #p.watson = 'A'
    p.dipole = jkext.convert.D2Cm([1.0, .0, .0])
    p.symmetry = 'N'
    for M in p.M:
        top = AsymmetricRotor(p, M, 0., 0.)
        for state in [State(0, 0, 0, 0, 0),
                      State(1, 0, 1, 0, 0), State(1, 1, 1, 0, 0), State(1, 1, 0, 0, 0),
                      State(2, 0, 2, 0, 0), State(2, 1, 2, 0, 0), State(2, 1, 1, 0, 0), State(2, 2, 1, 0, 0), State(2, 2, 0, 0, 0),
                      State(3, 0, 3, 0, 0), State(3, 1, 3, 0, 0), State(3, 1, 2, 0, 0), State(3, 2, 2, 0, 0),
                      State(3, 2, 1, 0, 0), State(3, 3, 1, 0, 0), State(3, 3, 0, 0, 0),]:
            print state.name(), "%10.3f %10.3g" % (jkext.convert.J2Hz(top.energy(state)) / 1e6, top.energy(state))
        top = AsymmetricRotor(p, M, 0., jkext.convert.kV_cm2V_m(100.))
#         for state in [State(0, 0, 0, 0, 0),
#                       State(1, 0, 1, 0, 0), State(1, 1, 1, 0, 0), State(1, 1, 0, 0, 0),
#                       State(2, 0, 2, 0, 0), State(2, 1, 2, 0, 0), State(2, 1, 1, 0, 0), State(2, 2, 1, 0, 0), State(2, 2, 0, 0, 0)]:
#             print state.name(), "%10.3f" % (jkext.convert.J2Hz(top.energy(state)) / 1e6,)
