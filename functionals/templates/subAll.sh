#!/bin/bash
source /nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/runpaths.source
{% if scan %}
source /nfs/slac/g/suncatfs/ksb/scan_env.bash
{% endif %}

/bin/hostname > host.txt

for d in strain*; do
  cd $d
  mpirun -n $LSB_MAX_NUM_PROCESSORS {{ exe }} &>run.log
  cd ..
