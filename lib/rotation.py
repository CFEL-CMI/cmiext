# -*- coding: utf-8; fill-column: 120 -*-
#
# This file is part of JK Python extensions
# Copyright (C) 2009 Jochen Küpper <software@jochen-kuepper.de>
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

"""Rotation routines"""

import scipy as num


def rotation_matrix(axis, angle):
    """Create 3D rotation matrix for rotation of |angle| around the space-fixed |axis| (x, y, z)."""
    from math import cos, sin
    from scipy import array
    if 'x' == axis:
        return array([[1, 0,           0],
                      [0, cos(angle), -sin(angle)],
                      [0, sin(angle),  cos(angle)]])
    elif 'y' == axis:
        return array([[cos(angle), 0, -sin(angle)],
                      [0,          1,  0],
                      [sin(angle), 0,  cos(angle)]])
    elif 'z' == axis:
        return array([[cos(angle), -sin(angle), 0],
                      [sin(angle),  cos(angle), 0],
                      [0,           0,          1]])
    else:
        raise "impossible axis specified"


def euler2cartesian(alpha, beta, gamma):
    """Conversion from Euler angles to cartesian coordinates -- not correct!

    Calculate the cartesian vector resulting from applying the specified Euler angles to the (0, 0, 1) unit vector.

    We chose to work in the ZXZ convension
    We use XYZ for the lab fixed system!!!
    this is from wikipedia they change between using xyz and XYZ for the fixed system
    o Rotate the XYZ-system about the z-axis by gamma. The X-axis is now at angle gamma with respect to the x-axis.
    o Rotate the XYZ-system again about the x-axis by beta. The Z-axis is now at angle beta with respect to the z-axis.
    o Rotate the XYZ-system a third time about the z-axis by alpha. The first and third axes are identical.
    this is equvelent to:
    (here the convension is changed)
    * Rotate the xyz-system about the z-axis by alpha. The x-axis now lies on the line of nodes.
    * Rotate the xyz-system again about the now rotated x-axis by beta. The z-axis is now in its final orientation,
      and the x-axis remains on the line of nodes.
    * Rotate the xyz-system a third time about the new z-axis by gamma.
    """
    from scipy import array, dot
    Z = num.array((0., 0., 1.), dtype=num.float64)
    R1 = rotation_matrix('z', alpha)
    R2 = rotation_matrix('x', beta)
    R3 = rotation_matrix('z', gamma)
    return num.dot(R3, num.dot(R2, num.dot(R1, Z)))


def euler2polar(alpha, beta, gamma):
    """Conversion from Euler angles to spherical polar coordinates.

    This converts uses euler2cartesian and then converts from cartesian to polar system.
    Returns r, theta, phi (angles in radians)"""
    from math import atan
    a = euler2cartesian(alpha, beta, gamma)
    r = (a[0]**2 + a[1]**2 + a[2]**2)**(1/2)
    theta = atan((a[0]**2+a[1]**2)**(1/2) / a[2])
    phi = atan(a[1] / a[0])
    return num.array((r, theta, phi))


def euler2polardegree(alpha, beta, gamma):
    """Conversion from Euler angles to spherical polar coordinates.

    Uses euler2polar, returns r, theta, phi (angles in degrees)"""
    return euler2polar(alpha, beta, gamma) * num.array((1, num.degrees(1), num.degrees(1)))


# some simple tests
if __name__ == "__main__":
    from math import pi
    print "\n\n"
    print euler2polardegree(0, pi/2, 0)
    print euler2polardegree(0, 0.1, 0.1)
    #
    print euler2cartesian(0, pi/2, 0)
    print euler2cartesian(pi/2, 0, 0)
    print euler2cartesian(0, pi, 0)
