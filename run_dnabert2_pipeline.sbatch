#!/bin/bash
#SBATCH --job-name=dnabert2_train
#SBATCH --output=logs/dnabert2_%j.out
#SBATCH --error=logs/dnabert2_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00

set -euo pipefail
bash ~/microbiome_taxonomy/run_dnabert2_pipeline.sh
