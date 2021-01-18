#!/bin/bash
source /nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/runpaths.source
source /nfs/slac/g/suncatfs/ksb/vasp54_env.bash

/bin/hostname > host.txt

for d in strain*/; do
  cd $d
  echo $d
  mpirun -n $LSB_MAX_NUM_PROCESSORS vasp &>rungga.log
  sed -i 's/GGA = PE/METAGGA = BF/g' INCAR
  mpirun -n $LSB_MAX_NUM_PROCESSORS vasp &>run.log
  rm -f WAVECAR
  rm -f CHG
  rm -f CHGCAR  cd ..
done
