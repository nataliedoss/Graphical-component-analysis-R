#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=sims
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --mail-user=natalie.doss@yale.edu

#module load Apps/R/3.5.1-generic
module load R
Rscript sims_plot.r 
