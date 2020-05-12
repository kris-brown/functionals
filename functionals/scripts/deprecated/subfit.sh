#!/bin/bash
source /nfs/slac/g/suncatfs/sw/py3.7.4/env.bash
source /nfs/slac/g/suncatfs/ksb/functionals/.env/bin/activate
export PYTHONPATH=/nfs/slac/g/suncatfs/ksb/functionals:$PYTHONPATH
python runfit.py &>run.log
