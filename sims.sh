#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=sims
#SBATCH --array=1-300
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mail-user=natalie.doss@yale.edu

module load R
Rscript sims.r ${SLURM_ARRAY_TASK_ID}
