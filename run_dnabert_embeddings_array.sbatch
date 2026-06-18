#!/usr/bin/env bash
#SBATCH --job-name=dnabert_embed
#SBATCH --partition=filemove
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/dnabert_embed_%A_%a.out
#SBATCH --error=logs/dnabert_embed_%A_%a.err

set -euo pipefail

source /apps/lib-osver/miniforge3/etc/profile.d/conda.sh
conda activate dnabert

IDS_FILE="${1:?Usage: sbatch --array=0-N run_dnabert_embeddings_array.sbatch ids.txt}"

ACC=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$IDS_FILE" | tr -d '[:space:]')

python ~/microbiome_taxonomy/extract_dnabert2_embeddings_from_fastq.py \
  --sample "$ACC" \
  --fastq_dir ~/scratch/fastq \
  --outdir ~/scratch/dnabert2_embeddings/"$ACC" \
  --batch_size 8 \
  --max_length 250
