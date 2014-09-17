#!/usr/bin/env tcsh
# install cmiext for CFEL-CMI use on BIRD (DESY compte cluster)
source /afs/desy.de/group/cfel/cfeld-cmi/bird/setup.csh
python setup.py install --home=${CMIBIRDPATH}
