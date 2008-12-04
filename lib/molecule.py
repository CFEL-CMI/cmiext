#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2008 Jochen Küpper <software@jochen-kuepper.de>
# see LICENSE file for details

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num
import numpy.linalg
import const
import tables

import starkeffect
from state import State


Masses = {'H': 1.0078250321, 'C': 12, 'N': 14.0030740052, 'O': 15.9949146221}
Ordernumbers = {'H': 1, 'D': 1, 'C': 6, 'N': 7, 'O': 8}


class Atom:
    """Representation of an atom

    Keeps a (private) list of Z and mass that can be accessed via public methods and a (public) copy of positions.

    Internal units are SI (i.e., m, kg, ...)! You can specify what length-unit is used on input (SI or Angstrom).
    """
    def __init__(self, symbol, position, length="SI"):
        self.__Z = Ordernumbers[symbol]
        self.__mass = Masses[symbol] * const.unified_atomic_mass
        self.__symbol = symbol
        self.position = num.array(position)
        assert(self.position.shape == (3,))
        if length == "Angstrom":
            self.position *= const.angstrom
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

    def __init__(self, atoms=None, storage=None, name="Generic molecule"):
        """Create Molecule from a list of atoms."""
        self.__atoms = atoms
        self.__name = name
        if storage != None:
            self.__storage = tables.openFile(storage, mode='a', title=name)
        else:
            self.__storage = None
        if atoms != None:
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
        rotcon = const.plancks_constant_h / (8 * num.pi**2 * eval[idx]) # in Hz
        # and provide corresponding eigenvectors of the axes (change columns!)
        axes = evec[:, idx]
        return rotcon, axes


    def rotate(self, rotation):
        """Rotate molecule using the right rotation matrix given (i.e., rotate all atomic coordinates)."""
        self.positions = num.dot(self.positions, rotation)
        return self.positions 


    def to_principal_axis_of_inertia(self):
        """Put molecule into its principal axis system."""
        # move to center of mass
        self.translate(-self.center_of_mass())
        # rotate to pricipal axis system
        rotcon, axes = self.principal_axis_of_inertia()
        self.rotate(axes)
        return rotcon


    def __get_starkeffect(self, state):
        entry = self.__storage.getNode("/dcstarkeffect/" + state.hdfname())
        data = entry.read(0)
        return data['field'][0], data['energy'][0]
    

    def __set_starkeffect(self, state, energies, fields):
        """Set potential |energies| for the specified |state| over the field-strength vector given in |fields|."""

        assert len(fields) == len(energies)

        class entry(tables.IsDescription):
            name = tables.StringCol(32)
            field = tables.Float64Col(len(fields))
            energy = tables.Float64Col(len(energies))
                
        try:
            group = self.__storage.root.dcstarkeffect
        except:
            group = self.__storage.createGroup('/', 'dcstarkeffect', 'DC Stark effect energy')
        try:
            table = self.__storage.removeNode("/dcstarkeffect/" + state.hdfname())
        except:
            pass
        table = self.__storage.createTable(group, state.hdfname(), entry, "potential energy")
        entry = table.row
        entry['name'] = state.name()
        entry['field'] = fields
        entry['energy'] = energies
        entry.append()
        self.__storage.flush()


    def starkeffect(self, state, energies=None, fields=None):
        """Get or set the potential energies as a function of the electric field strength.

        When |energies| and |fields| are None, return the Stark curve for the specified quantum state.

        When |energies| and |fields| are specified, save the Stark curve for the specified quantum state in the
        Molecule's HDF5 storage file.
        """
        if energies == None and fields == None:
            return self.__get_starkeffect(state)
        elif energies == None or fields == None:
            raise SyntaxError
        else:
            return self.__set_starkeffect(state, energies, fields)        


    def starkeffect_calculation(self, param):
        """Get all available energies from the given Starkeffect object and store them in our storage file."""
        energies = {}
        if 'A' == param.type:
            for M in param.M:
                for field in param.fields:
                    calc = starkeffect.AsymmetricRotor(param, M, field)
                    for state in calc.states():
                        id = state.id()
                        if energies.has_key(id):
                            energies[id].append(calc.energy(state))
                        else:
                            energies[id] = [calc.energy(state),]
        else:
            raise NotImplementedError("unknow rotor type in Stark energy calculation.")
        # store calculated values
        for id in energies.keys():
            self.starkeffect(State().fromid(id), energies[id], param.fields)


    def translate(self, translation):
        """Translate center of mass of molecule (i.e., translate all atoms)."""
        self.positions += translation
            


# some simple tests
if __name__ == "__main__":
    # test State
    state = State(0, 0, 0, 0)
    state.fromid(1000001)
    print state.name()
    # test Stark calculation and storage/retrieval
    from convert import *
    param = starkeffect.CalculationParameter
    param.rotcon = Hz2J(num.array([5e9, 2e9, 1.5e9]))
    param.dipole = D2Cm(num.array([1., 0., 0.]))
    mol = Molecule(storage="molecule.hdf")
    mol.starkeffect_calculation(param)
    for state in [State(0, 0, 0, 0), State(1, 0, 1, 1)]:
        fields, energies = mol.starkeffect(state)
        print V_m2kV_cm(fields), J2Hz(energies)
    


### Local Variables:
### fill-column: 132
### End:
