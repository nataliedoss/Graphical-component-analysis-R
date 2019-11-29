#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=test
#SBATCH --cpus-per-task=6
#SBATCH --time=24:00:00
#SBATCH --mail-user=natalie.doss@yale.edu

module load R
Rscript test.r
