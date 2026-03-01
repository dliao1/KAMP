#!/bin/bash
#SBATCH --array=1-24%4
#SBATCH --job-name=kamp_sim_large
#SBATCH --partition=week-long-cpu
#SBATCH --output=kamp_sim_large.out
#SBATCH --error=kamp_sim_large.err

module purge
module load R/4.3.2

export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3

JOBID=$SLURM_ARRAY_TASK_ID
Rscript kamp_sim_large.R $JOBID