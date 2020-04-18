CMI Python extensions
=====================

* License

The JK Python extensions are
  Copyright (C) 2008,2009,2016 Jochen Küpper
See License for details.


* Prerequisites

The extension package rests on the shoulders of giants, i.e., you need
a recent Python 3.x system and some other important Python packages,
including (but not necessarily limited to)
  - numpy
  - pytables

These requirements can all be met by using the Enthought Python
Distribution, for example.


* Installation

Adminstrator installation:
    python setup.py install

In order to install this extension module in user-space, run
    python setup.py install --user
and set your environment up for python to find it:
    setenv PYTHONUSERBASE=$HOME/.local

