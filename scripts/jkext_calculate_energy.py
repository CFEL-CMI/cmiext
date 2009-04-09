#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2009 Jochen Küpper <software@jochen-kuepper.de>
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
#
# some SGE commands for batch usage convenience
#$ -cwd
#$ -e $JOB_NAME.sge
#$ -o $JOB_NAME.sge
#$ -S $HOME/.python/bin/python
#$ -V
from __future__ import division

"""Calculate adiabatic energy curves and store them in HDF5 file

Copyright (C) 2008,2009 Jochen Küpper"""

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num
import getopt, sys

import jkext.molecule as molecule
import jkext as jk
import jkext.convert as convert
import jkext.molecule as molecule
import jkext.starkeffect as starkeffect
from jkext.state import State


def usage():
    # ToDo implement a useful usage description
    print "See script for details"


def main(args):
    try:
        opts, args = getopt.getopt(args[1:], "h", ["help", "Jmin=", "Jmax=", "dc-fields=",
                                                   "benzonitrile", "2,6-difluoro-iodobenzene"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    # default values
    Jmin = 0
    Jmax = 2
    param = starkeffect.CalculationParameter
    param.dcfields = convert.kV_cm2V_m(num.linspace(0., 100., 3))
    param.isomer = 0
    # scan commandline
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o == "--Jmin":
            Jmin = int(a)
        elif o == "--Jmax":
            Jmax = int(a)
        elif o == "--dc-fields":
            min, max, steps = a.split(":")
            param.dcfields = convert.kV_cm2V_m(num.linspace(float(min), float(max), int(steps)))
        elif o == "--benzonitrile":
            # Wohlfart, Schnell, Grabow, Küpper, J. Mol. Spec. 247, 119-121 (2008)
            name = "benzonitrile"
            param.watson = 'A'
            param.symmetry = 'C2a'
            param.rotcon = convert.Hz2J(num.array([5655.2654e6, 1546.875864e6, 1214.40399e6]))
            param.quartic = convert.Hz2J(num.array([45.6, 938.1, 500, 10.95, 628]))
            param.dipole = convert.D2Cm(num.array([4.5152, 0., 0.]))
        elif o == "--2,6-difluoro-iodobenzene":
            name = "2,6-difluoro-iodobenzene"
            param.watson = 'A'
            param.symmetry = 'C2a'
            param.rotcon = convert.Hz2J(num.array([1740e6, 713e6, 506e6]))
            param.quartic = convert.Hz2J(num.array([0., 0., 0., 0., 0.]))
            param.dipole = convert.D2Cm(num.array([2.25, 0., 0.]))
        else:
            assert False, "unhandled commandline option"

    # perform calculation
    mol = molecule.Molecule(storage=name+".hdf")
    mol.starkeffect_calculation(param)
    # close file
    del mol


if __name__ == "__main__":
    main(sys.argv)
