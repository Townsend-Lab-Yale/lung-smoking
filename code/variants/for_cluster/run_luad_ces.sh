#!/usr/bin/bash
#SBATCH --partition=pi_townsend
#SBATCH --job-name=luad_ces
#SBATCH --nodes=2
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=krishna.dasari@yale.edu

module load R

Rscript compute_selection.R