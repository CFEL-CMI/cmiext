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
# If you use this programm for scientific work, you must correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.
from __future__ import division

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

# really use scipy s numpy, so we are sure we use Fortran codes of eigvalsh and dgemm
import scipy as num
import scipy.linalg
import scipy.linalg.fblas
import jkext.convert
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
    fields = jkext.convert.kV_cm2V_m(num.linspace(0., 100., 2)) # V/m
    M = range(0, 2)
    Jmax_calc = 5
    Jmax_save = 2
    isomer = 0
    # molecular parameters
    rotcon = num.zeros((3,), num.float64)    # Joule
    quartic = num.zeros((5,), num.float64)   # Joule
    dipole = num.zeros((3,), num.float64)    # Coulomb meter
    watson=None
    symmetry='N'



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
        self.__stateorder_valid = False
        # save parameters internally
        self.__field = num.float64(field)
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
        # this requires a correct "global" value of self.__Jmin_matrixsize, which is set in __full_hamiltonian.
        # Therefore, we must be called only through __full_hamiltonian
        blockstart = J*(J-1) + J - self.__Jmin_matrixsize
        return blockstart + K + J


    def __recalculate(self):
        """Perform calculation of rotational state energies for current parameters"""
        self.__levels = {}
        blocks = self.__full_hamiltonian(self.__Jmin, self.__Jmax, self.__field, self.__symmetry)
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


    def __full_hamiltonian(self, Jmin, Jmax, field, symmetry):
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
        if None != field and self.__tiny < abs(field):
            self.__stark(hmat, Jmin, Jmax, field)
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


    def __stark(self, hmat, Jmin, Jmax, field):
        """Add the Stark-effect matrix element terms to hmat"""
        sqrt = num.sqrt
        M = self.__M
        muA, muB, muC = self.__dipole
        if self.__dipole_components[0]:
            # matrix elements involving µ_a
            for J in range(Jmin, Jmax):
                for K in range(-J, J+1):
                    if 0 != J:
                        hmat[self.__index(J, K), self.__index(J, K)] += -muA * field * M * K / (J*(J+1))
                    value = (-muA * field * sqrt((J+1)**2 - K**2) * sqrt((J+1)**2 - M**2)
                              / ((J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.__index(J+1, K), self.__index(J, K)] += value
                    hmat[self.__index(J, K), self.__index(J+1, K)] += value
            # final diagonal elements
            J = Jmax
            for K in range(-J, J+1):
                hmat[self.__index(J, K), self.__index(J, K)] += -1. * M * K / (J*(J+1)) * muA * field
        if self.__dipole_components[1]:
            # matrix elements involving µ_b
            for J in range(Jmin, Jmax):
                for K in range(-J, J+1):
                    if 0 != J:
                        value = -1 * M * muB * field * (sqrt((J-K) * (J+K+1) ) ) / (2*J*(J+1))
                        hmat[self.__index(J, K+1), self.__index(J, K)] += value
                        hmat[self.__index(J, K), self.__index(J, K+1)] += value
                    # J+1, K+1 / J-1, K-1 case
                    value = (muB * field * sqrt(((J+K+1) * (J+K+2)) * ((J+1)**2 - M**2))
                             / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.__index(J+1, K+1), self.__index(J, K)] += value
                    hmat[self.__index(J, K), self.__index(J+1, K+1)] += value
                    # J+1, K-1 / J-1, K+1 case
                    value = (-1 * muB * field * sqrt(((J-K+1) * (J-K+2)) * ((J+1)**2 - M**2))
                              / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.__index(J+1, K-1), self.__index(J, K)] += value
                    hmat[self.__index(J, K), self.__index(J+1, K-1)] += value
        if  self.__dipole_components[2]:
            # matrix elements involving µ_c
            for J in range(Jmin, Jmax):
                for K in range(-J, J+1):
                    if 0 != J:
                        value = 1j* M * muC * field * sqrt((J-K) * (J+K+1)) / (2*J*(J+1))
                        hmat[self.__index(J, K+1), self.__index(J, K)] += value
                        hmat[self.__index(J, K), self.__index(J, K+1)] += value
                    # J+1, K+1 / J-1, K-1 case
                    value = (-1j * muC * field * sqrt((J+K+1) * (J+K+2)) * sqrt((J+1)**2 - M**2)
                              / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.__index(J+1, K+1), self.__index(J, K)] += value
                    hmat[self.__index(J, K), self.__index(J+1, K+1)] += value
                    # J+1, K-1 / J-1, K+1 case
                    value = (-1j  * muC * field * sqrt((J-K+1) * (J-K+2)) * sqrt((J+1)**2 - M**2)
                              / (2*(J+1) * sqrt((2*J+1) * (2*J+3))))
                    hmat[self.__index(J+1, K-1), self.__index(J, K)] += value
                    hmat[self.__index(J, K), self.__index(J+1, K-1)] += value


    def __stateorder(self, symmetry):
        """Return a list with all states for the given |symmetry| and the current calculation parameters (Jmin, Jmax).

        Needs to be finished!
        """
        if False == self.__stateorder_valid:
            self.__stateorder_dict = {}
            M = self.__M
            iso = self.__isomer
            eigenvalues = {'Ep': [], 'Em': [], 'Op': [], 'Om': []}
            label = {'Ep': [], 'Em': [], 'Op': [], 'Om': []}
            for J in range(M, self.__Jmax+1):
                Ka = 0
                if 0 == J%2: # J even
                    for Kc in range(J,-1,-1):
                        if Ka%2 == 0 and Kc%2 == 0:  label['Ep'].append(State(J, Ka, Kc, M, iso))
                        elif Ka%2 == 0 and Kc%2 !=0: label['Em'].append(State(J, Ka, Kc, M, iso))
                        elif Ka%2 != 0 and Kc%2 ==0: label['Op'].append(State(J, Ka, Kc, M, iso))
                        else:                        label['Om'].append(State(J, Ka, Kc, M, iso)) # Ka%2 != 0 and Kc%2 != 0
                        if Kc > 0:
                            Ka = Ka+1
                            if Ka%2 == 0 and Kc%2 == 0:   label['Ep'].append(State(J, Ka, Kc, M, iso))
                            elif Ka%2 == 0 and Kc%2 != 0: label['Em'].append(State(J, Ka, Kc, M, iso))
                            elif Ka%2 != 0 and Kc%2 == 0: label['Op'].append(State(J, Ka, Kc, M, iso))
                            else:                         label['Om'].append(State(J, Ka, Kc, M, iso)) # Ka%2 != 0 and Kc%2 != 0
                else: # J odd
                    for Kc in range(J,-1,-1):
                        if Ka%2 ==0 and Kc%2 == 0:   label['Em'].append(State(J, Ka, Kc, M, iso))
                        elif Ka%2 ==0 and Kc%2 != 0: label['Ep'].append(State(J, Ka, Kc, M, iso))
                        elif Ka%2 !=0 and Kc%2 == 0: label['Om'].append(State(J, Ka, Kc, M, iso))
                        else:                        label['Op'].append(State(J, Ka, Kc, M, iso)) # Ka%2 != 0 and Kc%2 != 0
                        if Kc > 0:
                            Ka = Ka+1
                            if ((Ka%2==0) and  (Kc%2==0)):   label['Em'].append(State(J, Ka, Kc, M, iso))
                            elif ((Ka%2==0) and  (Kc%2!=0)): label['Ep'].append(State(J, Ka, Kc, M, iso))
                            elif ((Ka%2!=0) and  (Kc%2==0)): label['Om'].append(State(J, Ka, Kc, M, iso))
                            else:                            label['Op'].append(State(J, Ka, Kc, M, iso)) # Ka%2 != 0 and Kc%2 != 0
                # get block diagonal hamiltonian (make sure you calculate this in 'V'!)
                if 0 == J:
                    blocks = {'Ep': num.zeros((1, 1), num.float64)}
                else:
                    blocks = self.__full_hamiltonian(J, J, None, 'V')
                # store sorted eigenenergies for respective J and block
                for sym in blocks.keys():
                    if 0 < blocks[sym].size:
                        eigenvalues[sym] += num.sort(num.linalg.eigvalsh(num.array(blocks[sym]))).tolist()
            # sort assignments according to energy
            if 'V' == self.__symmetry:
                symmetries = ['Ep', 'Em', 'Op', 'Om']
            elif 'C2a' == self.__symmetry:
                symmetries = ['E', 'O']
                eigenvalues['E'] = eigenvalues['Ep'] + eigenvalues['Em']
                eigenvalues['O'] = eigenvalues['Op'] + eigenvalues['Om']
                label['E'] = label['Ep'] + label['Em']
                label['O'] = label['Op'] + label['Om']
            elif 'C2b' == self.__symmetry:
                raise NotImplementedError
            elif 'C2c' == self.__symmetry:
                raise NotImplementedError
            elif 'N' == self.__symmetry:
                symmetries = ['N']
                eigenvalues['N'] = eigenvalues['Ep'] + eigenvalues['Em'] + eigenvalues['Op'] + eigenvalues['Om']
                label['N'] = label['Ep'] + label['Em'] + label['Op'] + label['Om']
            else:
                raise NotImplementedError("Hamiltonian symmetry not implemented (yet)")
            for sym in symmetries:
                idx = num.argsort(eigenvalues[sym])
                self.__stateorder_dict[sym] = num.array(label[sym])[idx]
            self.__stateorder_valid = True
        return self.__stateorder_dict[symmetry]


    def __wang(self, hmat, symmetry, Jmin, Jmax):
        """Wang transform matrix and return a dictionary with the individual (sub)matrices."""
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
        # delete Wang matrix (it's not used anymore)
        del Wmat
        # sort out matrix blocks
        if 'N' == symmetry:
            # nothing to do, return
            blocks['N'] = hmat
        elif symmetry == 'V':
            # full Fourgroup symmetry
            # Only used for field free energies !
            # I^r representation, Wang transformed Hamiltonian factorizes into four submatrices E-, E+, O-, O+
            if Jmin != Jmax:
                raise NotImplementedError("Hamiltonian consists of more than one J-block. Is this field free ?")
            J = Jmax
            if J == 0:
                # return the single value in the transformed Hamiltonian
                blocks['Ep'] = num.array(hmat[0,0])
            else:
                blocksize = {}
                blocksize['Em'] = J//2
                blocksize['Ep'] = blocksize['Em'] + 1
                blocksize['Om'] = blocksize['Em'] + J%2
                blocksize['Op'] = blocksize['Om']
                blocks['Ep'] = num.zeros((blocksize['Ep'], blocksize['Ep']), num.float64)
                blocks['Em'] = num.zeros((blocksize['Em'], blocksize['Em']), num.float64)
                blocks['Op'] = num.zeros((blocksize['Op'], blocksize['Op']), num.float64)
                blocks['Om'] = num.zeros((blocksize['Om'], blocksize['Om']), num.float64)
                #fill E-
                if blocksize['Em'] > 0: # for J == 1, E- does not exist
                    for i in range(blocksize['Em']):
                        column = 2*i + J%2
                        for m in range(blocksize['Em']):
                            row = 2*m + J%2
                            blocks['Em'][m, i] = hmat[row, column]
                # fill E+
                for i in range(blocksize['Ep']):
                    column = 2*i + J
                    for m in range(blocksize['Ep']):
                        row = 2*m + J
                        blocks['Ep'][m, i] = hmat[row, column]
                # fill O-
                for i in range(blocksize['Om']):
                    column = 2*i - J%2 + 1
                    for m in range(blocksize['Om']):
                        row = 2*m - J%2 + 1
                        blocks['Om'][m, i] = hmat[row, column]
                # fill O+
                for i in range(blocksize['Op']):
                    column = 2*i + J + 1
                    for m in range(blocksize['Op']):
                        row = 2*m + J + 1
                        blocks['Op'][m, i] = hmat[row, column]
        elif symmetry == 'C2a':
            # C2 rotation about a-axis is symmetry element
            #
            # I^r representation, Wang transformed Hamiltonian factorizes into two submatrices E (contains E- and
            # E+) and O (contains O- and O+). In this case E and O corresponds to columns with Ka even and odd,
            # respectively.
            matrixsize_Jmin = Jmin *(Jmin-1) + Jmin
            matrixsize = (Jmax + 1) * Jmax + Jmax + 1 - matrixsize_Jmin
            matsize = matrixsize
            ESize = 0
            OSize = 0
            if (Jmax-Jmin)%2 != 0: #even number of J-blocks
                ESize = matsize/2
                OSize = matsize/2
            elif ((Jmin%2 != 0) and  (Jmax%2 !=0)): # odd number of J-blocks, Jmin and Jmax are both odd
                ESize = (matsize-1)/2
                OSize = (matsize+1)/2
            else: #odd number of J-blocks, Jmin and Jmax are both even
                ESize = (matsize+1)/2
                OSize = (matsize-1)/2
            Even = num.zeros((ESize, ESize), num.float64)
            Odd = num.zeros((OSize, OSize), num.float64)
            if 0 == Jmin%2: # start with even J block
                for i in range(0, matsize, 2):
                    for m in range(0,matsize,2):
                        Even[i/2,m/2] = hmat[i,m]
                for i in range(1, matsize, 2):
                    for m in range(1, matsize, 2):
                        Odd[(i-1)/2, (m-1)/2] = hmat[i,m]
            else: # start with odd J block
                for i in range(0, matsize, 2):
                    for m in range(0,matsize,2):
                        Odd[i/2,m/2] = hmat[i,m]
                for i in range(1, matsize, 2):
                    for m in range(1, matsize, 2):
                        Even[(i-1)/2, (m-1)/2] = hmat[i,m]
            blocks['E'] = Even
            blocks['O'] = Odd
        elif symmetry == 'C2b': # C2 rotation about b-axis is symmetry element
            raise NotImplementedError("Hamiltonian symmetry 'b' not implemented yet")
        elif symmetry == 'C2c': # C2 rotation about c-axis is symmetry element
            raise NotImplementedError("Hamiltonian symmetry 'c' not implemented yet")
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
    p.Jmax_calc = 2
    p.rotcon = jkext.convert.Hz2J(num.array([5e9, 2e9, 1.5e9]))
    p.quartic = jkext.convert.Hz2J([1e3, 1e3, 1e3, 1e3, 1e3])
    p.watson = 'A'
    p.dipole = jkext.convert.D2Cm([1., 1., 1.])
    top = AsymmetricRotor(p, 0, 0)
    for state in [State(0, 0, 0, 0, 0),
                  State(1, 0, 1, 0, 0), State(1, 1, 1, 0, 0), State(1, 1, 0, 0, 0),
                  State(2, 1, 2, 0, 0)]:
        print state.name(), "%10.3f" % (jkext.convert.J2Hz(top.energy(state)) / 1e6,)
    top = AsymmetricRotor(p, 0, jkext.convert.kV_cm2V_m(100.))
    for state in [State(0, 0, 0, 0, 0),
                  State(1, 0, 1, 0, 0), State(1, 1, 1, 0, 0), State(1, 1, 0, 0, 0),
                  State(2, 1, 2, 0, 0)]:
        print state.name(), "%10.3f" % (jkext.convert.J2Hz(top.energy(state)) / 1e6,)
