#!/bin/bash
#SBATCH --array=1-24%4
#SBATCH --job-name=kamp_expec_job
#SBATCH --partition=day-long-cpu
#SBATCH --output=kamp_expec_job.out
#SBATCH --error=kamp_expec_job.err

module purge
module load R/4.3.2

export R_LIBS_USER=$HOME/R/4.3.2

JOBID=$SLURM_ARRAY_TASK_ID
Rscript sim_expectation.R $JOBID