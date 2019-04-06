#!/bin/bash
#source /nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/runpaths.source

/bin/hostname > host.txt

mpirun -n $LSB_MAX_NUM_PROCESSORS /nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/vasp.5.4.1/bin/vasp_gam &> run.log

mv INCAR2 INCAR

mpirun -n $LSB_MAX_NUM_PROCESSORS /nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/vasp.5.4.1/bin/vasp_gam &>> run.log

rm WAVECAR
