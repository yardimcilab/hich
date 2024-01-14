#!/bin/bash
#SBATCH --partition exacloud
#SBATCH --job-name HiCkory
#SBATCH -c 24
#SBATCH --mem 25G
#SBATCH --output=logs/sbatch/HiCkory_%j.out         ### File in which to store job output
#SBATCH --error=logs/sbatch/HiCkory_%j.err          ### File in which to store job error messages
#SBATCH --time=10:00:00       ### Wall clock time limit in Days-HH:MM:SS
mkdir -p logs/sbatch
source setup_env.sh
source snakemake.sh
