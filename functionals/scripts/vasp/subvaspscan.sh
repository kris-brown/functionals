#!/bin/tcsh
source /nfs/slac/g/suncatfs/sw/vasp/vbfscan/setupenv
echo ${{HOST}} > host.txt

mpirun -n ${{LSB_MAX_NUM_PROCESSORS}} vasp > run1.log

cp OUTCAR OUTCAR0
cp INCAR INCAR0
cp INCAR2 INCAR

mpirun -n ${{LSB_MAX_NUM_PROCESSORS}} vasp > run2.log

rm -f WAVECAR
