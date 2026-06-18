#!/usr/bin/env bash
#SBATCH --job-name=sra_host_meta_kraken
#SBATCH --partition filemove
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_%x_%A_%a.out
#SBATCH --error=slurm_%x_%A_%a.err
#SBATCH --array=0-0

# IMPORTANT:
# Update the array range after you know how many IDs you have (see below)

set -euo pipefail

IDS_FILE="${1:?Usage: sbatch submit.sbatch ids.txt}"

# If your cluster uses modules, you might need this here instead of in the worker:
# module purge
# module load miniconda

# Run one sample per array task
bash run_one_sample.sh "$IDS_FILE"

