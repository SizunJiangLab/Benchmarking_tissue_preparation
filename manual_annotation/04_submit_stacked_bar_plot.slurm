#!/bin/bash
#SBATCH --job-name=figure1
#SBATCH --partition=batch_cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=output_%A_%a.log
#SBATCH --error=error_%A_%a.log
#SBATCH --time=03:00:00

ml Miniforge3
conda activate clustering_3.11

python ./04_stacked_bar_plots.py
