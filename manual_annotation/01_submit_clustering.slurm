#!/bin/bash
#SBATCH --job-name=clustering
#SBATCH --partition=batch_gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=output_%A_%a.log
#SBATCH --error=error_%A_%a.log
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --array=1-24

ml Miniforge3
conda activate clustering_3.11

python ./01_clustering.py ${SLURM_ARRAY_TASK_ID}
