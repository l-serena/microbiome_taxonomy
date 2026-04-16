#!/bin/bash
#SBATCH --job-name=dnabert2_train
#SBATCH --output=logs/dnabert2_%j.out
#SBATCH --error=logs/dnabert2_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

set -euo pipefail

echo "=== Job started on $(hostname) ==="
date

# --- conda init ---
source /apps/lib-osver/miniforge3/etc/profile.d/conda.sh
conda activate dnabert

echo "=== Python ==="
which python
python --version

echo "=== STEP 1: build dataset ==="
python ~/microbiome_taxonomy/build_dnabert2_read_dataset_split.py

echo "=== STEP 2: train DNABERT2 ==="
python ~/microbiome_taxonomy/train_dnabert2_read_classifier_split.py \
  --train_csv ~/scratch/dnabert2_reads_10samples.train.csv \
  --test_csv ~/scratch/dnabert2_reads_10samples.test.csv \
  --output_dir ~/scratch/dnabert2_read_model_10samples \
  --epochs 3 \
  --batch_size 8 \
  --min_class_count 20

echo "=== DONE ==="
date
