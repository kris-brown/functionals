#!/bin/bash
#source /nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/runpaths.source

/bin/hostname > host.txt

mpirun -n $LSB_MAX_NUM_PROCESSORS /nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/vasp.5.4.1beefcar/bin/vasp_std &>run.log

mv INCAR INCAR0
mv INCAR2 INCAR

mpirun -n $LSB_MAX_NUM_PROCESSORS /nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/vasp.5.4.1beefcar/bin/vasp_std &>run.log

rm -f WAVECAR
