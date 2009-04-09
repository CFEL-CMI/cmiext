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

"""Plot adiabatic energy curves from HDF5 files to matplotlib graph

Copyright (C) 2009 Jochen Küpper"""

__author__ = "Jochen Küpper <software@jochen-kuepper.de>"

import numpy as num
import getopt, sys
import matplotlib.pyplot as plt

import jkext.molecule as molecule
from jkext.convert import *
from jkext.state import State


def usage():
    # ToDo implement a useful usage description
    print "See script for details"


def main(args):
    try:
        opts, args = getopt.getopt(args[1:], "h", ["help", "Jmin=", "Jmax="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    # default values
    Jmin = 0
    Jmax = 2
    # scan commandline
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-s", "--storage"):
            storagename = a
        elif o == "--Jmin":
            Jmin = int(a)
        elif o == "--Jmax":
            Jmax = int(a)
        else:
            assert False, "unhandled commandline option"

    plt.clf()
    # loop over all remaining arguments -- asumming its filenames of HDF5 Stark-files
    for name in args:
        print name
        # create Molecule object and specify storage file
        mol = molecule.Molecule(storage=name)

        for J in range(Jmin, Jmax+1):
            Ka = 0
            for Kc in range(J, -1, -1):
                state = State(J, Ka, Kc, 0, 0)
                fields, energies = mol.starkeffect(state)
                plt.plot(V_m2kV_cm(fields), J2Hz(energies) / 1e6)
                if Kc > 0:
                    Ka += 1
                    state = State(J, Ka, Kc, 0, 0)
                    fields, energies = mol.starkeffect(state)
                    plt.plot(V_m2kV_cm(fields), J2Hz(energies) / 1e6)

    plt.show()


if __name__ == "__main__":
    main(sys.argv)
