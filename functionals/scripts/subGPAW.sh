#!/bin/bash
#SBATCH -p suncat,iric
#################
#SBATCH --job-name=myjob
#################
#a file for job output, you can check job progress
#SBATCH --output=myjob.out
#################
# a file for errors from the job
#SBATCH --error=myjob.err
#################
#SBATCH --time={Job.time:d}:55:00
#################
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-core=1
#################
#SBATCH --mem-per-cpu=4000
#################
#SBATCH --mail-type=FAIL
#################
#SBATCH  --mail-user=ksb@stanford.edu
#################
#figure out how many cores have been requested with the above
#slurm parameters

#load gpaw-specific paths
source /scratch/users/ksb/gpaw/paths.bash
#run parallel gpaw
mpirun -n $SLURM_TASKS_PER_NODE gpaw-python gpawrun.py
