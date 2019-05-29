#!/bin/bash
echo $HOST > host.txt

mpirun -n $LSB_MAX_NUM_PROCESSORS /nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/vasp.5.4.1/bin/vasp_std &> run.log

mv INCAR INCAR0
mv INCAR2 INCAR

vasp-ver-bsub-native bfscan -n 16 -o qlog.txt -q suncat3 -W 20:00

rm -f WAVECAR
