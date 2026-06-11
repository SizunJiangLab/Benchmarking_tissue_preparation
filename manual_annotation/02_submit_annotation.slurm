#!/bin/bash
#SBATCH --job-name=annotation
#SBATCH --partition=batch_cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=output_%A_%a.log
#SBATCH --error=error_%A_%a.log
#SBATCH --time=03:00:00

ml Miniforge3
conda activate clustering_3.11

# Data root for input/, output/, and cluster_annotations.json.
# Override by exporting REGISTERED_REPORT_DIR before sbatch; defaults to the cluster path.
export REGISTERED_REPORT_DIR="${REGISTERED_REPORT_DIR:-/registered_report}"

python ./02_annotation.py
