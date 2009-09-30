/* -*- coding: utf-8; fill-column: 120 -*-

Copyright (C) 2009 Jochen Küpper <python@jochen-kuepper.de>

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <Python.h>
#include "wigner.h"

static PyObject *_wigner_drot(PyObject *self, PyObject *args)
{
    int j, m1, m2;
    double theta, result;
    if (!PyArg_ParseTuple(args, "iiid", &j, &m1, &m2, &theta))
        return NULL;
    result = gsl_sf_wigner_drot(j, m1, m2, theta);
    return Py_BuildValue("d", result);
}


static PyMethodDef _wigner_methods[] = {
    {"drot", _wigner_drot, METH_VARARGS, "Reduced Wigner d rotation matrix."},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC init_wigner(void)
{
    (void) Py_InitModule("_wigner", _wigner_methods);
}
