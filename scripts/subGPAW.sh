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
#SBATCH --time={Job.time}:55:00
#################
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#################
#SBATCH --mem-per-cpu=4000
#################
#SBATCH --mail-type=FAIL
#################
#SBATCH  --mail-user=ksb@stanford.edu
#################
#figure out how many cores have been requested with the above
#slurm parameters

NTASKS=`echo $SLURM_TASKS_PER_NODE|tr '(' ' '|awk '{{print $1}}'`
NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`
NCPU=`echo " $NTASKS * $NNODES " | bc`
echo "NTASKS = $NTASKS"
echo "NNODES = $NNODES"
echo "NCPU = $NCPU"
#load gpaw-specific paths
source /scratch/users/ksb/gpaw/paths.bash
#run parallel gpaw
mpirun -n $NCPU gpaw-python gpawrun.py
