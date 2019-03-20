#!/bin/bash
source /nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/runpaths.source

echo $HOST > host.txt

mpirun -n $LSB_MAX_NUM_PROCESSORS /nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/vasp.5.4.1/bin/vasp_std &>run.log
