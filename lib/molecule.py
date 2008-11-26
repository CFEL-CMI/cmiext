#!/usr/bin/env python
#
# Copyright (C) 2008 Jochen Küpper

import numpy as num
import numpy.linalg
import pygsl.const

Masses = {'H': 1.0078250321, 'C': 12, 'N': 14.0030740052, 'O': 15.9949146221}
Ordernumbers = {'H': 1, 'D': 1, 'C': 6, 'N': 7, 'O': 8}


class Atom:
    """Representation of an atom

    Keeps a (private) list of Z and mass that can be accessed via public methods and a (public) copy of positions.

    Internal units are SI (i.e., m, kg, ...)! You can specify what length-unit is used on input (SI or Angstrom).
    """
    def __init__(self, symbol, position, length="SI"):
        self.__Z = Ordernumbers[symbol]
        self.__mass = Masses[symbol] * pygsl.const.unified_atomic_mass
        self.__symbol = symbol
        self.position = num.array(position)
        assert(self.position.shape == (3,))
        if length == "Angstrom":
            self.position *= pygsl.const.angstrom
        else:
            assert(length == "SI")

    def mass(self):
        return self.__mass

    def ordernumber(self):
        return self.__Z

    def position(self):
        return self.position

    def symbol(self):
        return self.__symbol

    def Z(self):
        return self.__Z



class Molecule:
    """Representation of a Molecule

    Keeps a (private) list of atoms and masses that can be accessed via |atom| and |masses| and a (public) array of
    positions.
    """

    def __init__(self, atoms=[]):
        """Create Molecule from a list of atoms."""
        self.__atoms = atoms
        self.__update()

    def __update(self):
        self.masses = num.zeros((len(self.__atoms),))
        for i in range(self.masses.shape[0]):
            self.masses[i] = self.__atoms[i].mass()
        self.positions = num.zeros((len(self.__atoms), 3))
        for i in range(self.positions.shape[0]):
            self.positions[i,:] = self.__atoms[i].position


    def center_of_mass(self):
        """Calculate center of mass of molecule"""
        return num.sum(num.outer(self.masses, [1,1,1]) * self.positions, axis=0) / num.sum(self.masses)


    def plot(self):
        """Create 2D plot of molecule."""
        import matplotlib.pyplot as plt
        plt.plot(self.positions[:,0], self.positions[:,1], 'o')
        plt.show()


    def principal_axis_of_inertia(self):
        """Calulate principal axes of inertia and the corresponding rotational constants.

        Returns array of rotational constants [A, B, C] (Hz) and the right eigenvectors of the inertial tensor.
        """
        inertial_tensor = num.zeros((3,3))
        for i in range(0,3):
            for j in range(0,3):
                for k in range(len(self.masses)):
                    m = self.masses[k]
                    r = self.positions[k]
                    if i == j:
                        r2 = num.dot(r, r) # square of length
                    else:
                        r2 = 0.
                    inertial_tensor[i,j] += m * (r2 - r[i]*r[j])
        eval, evec = num.linalg.eigh(inertial_tensor)
        # sort eigenvalues and eigenvetors
        idx = num.argsort(eval) # sort moments of inertia in ascending order
        # calculate rotational constants in Hz
        rotcon = pygsl.const.plancks_constant_h / (8 * num.pi**2 * eval[idx]) # in Hz
        # and provide corresponding eigenvectors of the axes (change columns!)
        axes = evec[:, idx]
        return rotcon, axes


    def rotate(self, rotation):
        """Rotate molecule using the right rotation matrix given (i.e., rotate all atomic coordinates)."""
        num.dot(self.positions, rotation)


    def to_principal_axis_of_inertia(self):
        """Put molecule into its principal axis system."""
        # move to center of mass
        self.translate(-self.center_of_mass())
        # rotate to pricipal axis system
        rotcon, axes = self.principal_axis_of_inertia()
        self.rotate(axes)
        return rotcon
            

    def translate(self, translation):
        """Translate center of mass of molecule (i.e., translate all atoms)."""
        self.positions += translation
            


    


### Local Variables:
### fill-column: 120
### End:
