#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2008 Jochen Küpper <software@jochen-kuepper.de>
# see LICENSE file for details

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num
import numpy.linalg
from molecule import State


class CalculationParameter:
    """Container of parameters for calculation of Stark energies.

    Calculate energy for the specified |fields| (V/m) and rotor type:
    - 'A': asymmetric top

    The following parameter are used for
    - an asymmetric top: M, Jmin, Jmax_calc, Jmax_save
    """
    type = 'A'
    fields = num.linspace(0., 100., 5) # V/m
    M = range(0, 2)
    Jmin = 0
    Jmax_calc = 5
    Jmax_save = 2
    isomer = 0
    # molecular parameters
    rotcon = num.zeros((3,))            # MHz
    quartic = num.zeros((5,))           # MHz
    dipole = num.zeros((3,))            # D
    reduction=None
    symmetry=None



class AsymmetricRotor:
    """Representation of an asymmetric top for energy level calculation purposes.

    This object will caclulate rotational energies at the specified field strength for the given M-value and J-range and
    all K's.
    """

    def __init__(self, param, M, field):
        """Save the relevant parameters"""
        assert 'A' == param.type
        # we have not yet calculated the correct energies - mark invalid
        self.__valid = False
        # save parameters internally
        self.__field = float(field)
        self.__rotcon = num.array(param.rotcon)
        self.__quartic = num.array(param.quartic)
        self.__dipole = num.array(param.dipole)
        self.__reduction = param.reduction
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
        for J in range(self.__Jmin, self.__Jmax+1):
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
        self.__levels = {}
        self.__full_hamiltonian()
        blocks = self.__wang()
        for symmetry in blocks.keys():
            eval = num.linalg.eigvalsh(blocks[symmetry]) # calculate only energies
            eval = num.sort(eval)
            if self.__symmetry == None:
                i = 0
                for stateid in self.__stateidorder(symmetry):
                    self.__levels[stateid] = eval[i]
                    i += 1
            else:
                raise NotImplementedError("Symmetric cases not implemented yet")
        # done - data is now valid
        self.__valid == True
        

    def __full_hamiltonian(self):
        # create hamiltonian matrix
        self.__hmat = num.zeros((self.__matrixsize, self.__matrixsize))
        # start matrix with appropriate field-free rigid-rotor terms
        self.__rigid()
        # add appropriate field-free centrifugal distortion terms
        if self.__reduction == 'A':
            raise NotImplementedError("Watson's A-reduction is not implemented (yet)")
            self.__watson_A()
        elif self.__reduction == 'S':
            raise NotImplementedError("Watson's S-reduction is not implemented (yet)")
            self.__watson_S()
        else:
            assert self.__reduction == None
        # fill matrix with appropriate Stark terms
        self.__stark()

        
    def __rigid(self):
        """Add the rigid-rotor matrix element terms to self.__hmat"""
        sqrt = num.sqrt
        A, B, C = self.__rotcon.tolist()
        for J in range(self.__Jmin, self.__Jmax+1):
            for K in range(-J, J+1):
                self.__hmat[self.__index(J, K), self.__index(J, K)] = (B+C)/2 * (J*(J+1) - K**2) + A * K**2
            for K in range (-J, J-2+1):
                self.__hmat[self.__index(J, K+2), self.__index(J, K)] \
                    = (C-B)/4 * sqrt(J*(J+1) - (K+1)*(K+2)) * sqrt(J*(J+1) - K*(K+1))
                self.__hmat[self.__index(J, K), self.__index(J, K+2)] \
                    = (C-B)/4 * sqrt(J*(J+1) - (K+1)*(K+2)) * sqrt(J*(J+1) - K*(K+1))


    def __stark(self):
        """Add the Stark-effect matrix element terms to self.__hmat."""
        pass



    def __stateidorder(self, symmetry):
        """Return a list with all stateids for the given |symmetry| and the current calculation parameters 
        (Jmin, Jmax_save).""" 
        ids = []
        for state in self.states():
            ids.append(state.id())
        return ids        
        

    def __wang(self):
        """Wang transform matrix and return a dictionary with the individual submatrices."""
        blocks = {}
        if self.__symmetry == None:
            # nothing to do, return
            blocks['N'] = self.__hmat
        else:
            # set up Wang matrix
            Wmat = num.zeros(self.__hmat.shape)
            value = 1/num.sqrt(2.)
            for J in range(self.__Jmin, self.__Jmax + 1):
                for K in range(-J, 0):
                    Wmat[self.__index(J,  K), self.__index(J,  K)] = -value
                    Wmat[self.__index(J, -K), self.__index(J,  K)] = value
                    Wmat[self.__index(J,  K), self.__index(J, -K)] = value
                    Wmat[self.__index(J, -K), self.__index(J, -K)] = value
                Wmat[self.__index(J, 0), self.__index(J, 0)] = 1.
            # transform Hamiltonian matrix
            print self.__hmat
            self.__hmat = num.dot(num.dot(Wmat, self.__hmat), Wmat)
            print self.__hmat
            # sort out matrix blocks
            if self.__symmetry == 'V': 
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
        pass


    def __watson_S(self):
        """Add the centrifugal distortion matrix element terms in Watson's S reduction to self.__hmat."""
        pass



# some simple tests
if __name__ == "__main__":
    print
    p = CalculationParameter
    p.rotcon = num.array([5e9, 2e9, 1.8e9])
    top = AsymmetricRotor(p, 0, 0)
    print top.energy(State(1, 0, 1, 0, 0)) / 1e6


### Local Variables:
### fill-column: 120
### End:
