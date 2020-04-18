# CMI Python extensions

## License

The CMI Python extensions are _Copyright (C) 2008,2009,2016,2020 Jochen Küpper_; see LICENSE for details.


## Prerequisites

The extension package rests on the shoulders of giants, i.e., you need
a recent Python 3.x system and some other important Python packages,
including (but not necessarily limited to)
  - numpy
  - pytables


## Installation

Adminstrator installation:
```shell
python setup.py install
```

In order to install this extension module in user-space, [set up your
environment](https://docs.python.org/3/using/cmdline.html#envvar-PYTHONUSERBASE)
for python to find it, e.g.,
```shell
setenv PYTHONUSERBASE=$HOME/.local
```
and run the install command
```shell
python setup.py install --user
```


<!-- Put Emacs local variables into HTML comment
Local Variables:
coding: utf-8
fill-column: 80
End:
-->
