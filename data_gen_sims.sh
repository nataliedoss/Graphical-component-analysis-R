#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=data_gen_sims
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=natalie.doss@yale.edu

module load R
Rscript data_gen_sims.r 
