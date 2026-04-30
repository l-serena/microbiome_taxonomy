#!/bin/bash
#SBATCH --job-name=dnabert2_200k
#SBATCH --output=/PHShome/sl1172/microbiome_taxonomy/logs/dnabert2_200k_%j.out
#SBATCH --error=/PHShome/sl1172/microbiome_taxonomy/logs/dnabert2_200k_%j.err
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

set -euo pipefail

source /apps/lib-osver/miniforge3/etc/profile.d/conda.sh
conda activate dnabert

echo "=== HOST ==="
hostname
date

echo "=== PYTHON ==="
which python
python --version

echo "=== TRAINING ==="
python ~/microbiome_taxonomy/train_dnabert2_read_classifier_split.py \
  --train_csv ~/scratch/dnabert2_reads_10samples.sampled_200k.train.csv \
  --test_csv ~/scratch/dnabert2_reads_10samples.sampled_200k.test.csv \
  --output_dir ~/scratch/dnabert2_read_model_200k \
  --epochs 1 \
  --batch_size 8 \
  --min_class_count 20 \
  --eval_steps 1000 \
  --save_steps 1000

echo "=== DONE ==="
date
