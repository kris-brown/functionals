#!/bin/bash
#SBATCH -p suncat
#################
#SBATCH --job-name=myjob
#################
#a file for job output, you can check job progress
#SBATCH --output=myjob.out
#################
# a file for errors from the job
#SBATCH --error=myjob.err
#################
#SBATCH --time={Job.time:d}:03:00
#################
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-core=1
#################
#SBATCH --mem=50000
#################
#SBATCH --mail-type=FAIL
#################
#SBATCH  --mail-user=ksb@stanford.edu
#################
#figure out how many cores have been requested with the above
#slurm parameters


source /home/users/ksb/vossj/vasp/nusherlock/source.txt
export OMP_NUM_THREADS=1
#run parallel vasp
mpirun -n $SLURM_TASKS_PER_NODE /home/users/ksb/vossj/vasp/nusherlock/vasp.5.4.1/bin/vasp_std
