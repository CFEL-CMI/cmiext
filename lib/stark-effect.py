#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2008 Jochen Küpper

import numpy as num
import numpy.linalg
import pygsl.const


class State:
    """State label"""

    def __init__(self, J, Ka, Kc, M, isomer=0):
        self.labels = num.array([J, Ka, Kc, M, isomer])

    def toarray(self):
        return self.labels

    def tolist(self):
        return self.labels.tolist()


class AsymmetricRotor:
    """Representation of an asymmetric top for energy level calculation purposes."""

    def __init__(self, rotcon, quartic=[0., 0., 0., 0., 0.], dipole=None, reduction=None, symmetry=None):
        self.__valid = False # we have not yet calculated the correct energies
        self.__rotcon = num.array(rotcon)
        assert self.__rotcon.shape == (3,)
        self.__quartic = num.array(quartic)
        assert self.__quartic.shape == (5,)
        self.__reduction = reduction
        self.__symmetry = symmetry

        
    def energy(self, state):
        """Return Stark energy for |state|."""
        return self.level(state)[0]
           

    def level(self, state):
        """Return Stark energy and eigenvector for |state|."""
        if not self.__valid:
            self.__recalculate()
        return self.__levels[state]
           

    def quantum_numbers(self, M=0, Jmin=0, Jmax=10, Jmax_save=2):
        self.__M = M
        self.__Jmin = Jmin
        self.__Jmax = Jmax
        self.__Jmax_save = Jmax_save
        self.__matrixsize_Jmin = self.__Jmin *(self.__Jmin-1) + self.__Jmin
        self.__matrixsize = (self.__Jmax + 1) * self.__Jmax + self.__Jmax + 1 - self.__matrixsize_Jmin


    def __index(self, J, K):
        blockstart = J*(J-1) + J - self.__matrixsize_Jmin
        return blockstart + K + J
        

    def __recalculate(self):
        self.__full_hamiltonian()
        print self.__hmat
        eval, evec = num.linalg.eigh(self.__hmat)
        idx = num.argsort(eval)
        print eval
        print evec
        eval, evec = eval[idx], evec[:, idx]
        self.__levels = {}
        self.__levels[State(0, 0, 0, 0, 0)] = [eval[0], evec[:,0]]
        self.__levels[State(1, 0, 1, 0, 0)] = [eval[1], evec[:,1]]
        print self.__levels
        
    def __full_hamiltonian(self):
        # create hamiltonian matrix
        self.__hmat = num.zeros((self.__matrixsize, self.__matrixsize))
        # fill matrix with appropriate field-free terms
        self.__rigid()
        if self.__reduction == 'A':
            self.__watson_A()
        elif self.__reduction == 'S':
            print "not implemented"
            exit(1)
        else:
            assert self.__reduction == None
        # fill matrix with appropriate Stark terms
        self.__stark()

        
    def __rigid(self):
        """Add the rigid-rotor matrix element terms to self.__hmat."""
        from math import sqrt
        A, B, C = self.__rotcon.tolist()
        print A, B, C
        for J in range(self.__Jmin, self.__Jmax):
            for K in range(-J, J):
                self.__hmat[self.__index(J, K), self.__index(J, K)] = (B+C)/2 * (J*(J+1) - K**2) + A * K**2
            for K in range (-J, J-2):
                self.__hmat[self.__index(J, K+2), self.__index(J, K)] = (C-B)/4 * sqrt(J*(J+1) - (K+1)*(K+2)) * sqrt(J*(J+1) - K*(K+1))
                self.__hmat[self.__index(J, K), self.__index(J, K+2)] = (C-B)/4 * sqrt(J*(J+1) - (K+1)*(K+2)) * sqrt(J*(J+1) - K*(K+1))


    def __stark(self):
        pass

    def __watson_A(self):
        pass




if __name__ == "__main__":
    # test calculation for benzonitrile (J Mol Spec ...)
    rotcon = num.array([5e9, 2e9, 1.8e9])
    top = AsymmetricRotor(rotcon)
    top.quantum_numbers(0, 0, 1, 1)
    state = State(0, 0, 0, 0, 0)
    print top.energy(state)


### Local Variables:
### fill-column: 120
### End:
