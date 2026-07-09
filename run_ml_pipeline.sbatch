#!/usr/bin/env bash
#SBATCH --job-name=ml_taxonomy
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/ml_taxonomy_%j.out
#SBATCH --error=logs/ml_taxonomy_%j.err

set -euo pipefail

REPO="${REPO:-$HOME/microbiome_taxonomy}"
source /apps/software/Miniforge3/24.11.3-0/etc/profile.d/conda.sh

mkdir -p "$REPO/logs" ~/scratch/models

echo "=== Pipeline status ==="
bash "$REPO/check_pipeline_status.sh" "$REPO/ids.txt"

echo "=== Build merged dataset ==="
conda activate dnabert
python "$REPO/build_merged_taxonomy_dataset.py" \
  --ids_file "$REPO/ids.txt" \
  --max_reads_per_sample 5000

echo "=== Train hierarchical greedy classifier ==="
python "$REPO/train_hierarchical_classifier.py" \
  --datadir ~/scratch/merged_taxonomy_ml \
  --outdir ~/scratch/models/hierarchical_greedy \
  --min_class_count 20

echo "=== Train read GNN (species) ==="
python "$REPO/train_read_gnn.py" \
  --datadir ~/scratch/merged_taxonomy_ml \
  --outdir ~/scratch/models/read_gnn_species \
  --label_rank species \
  --epochs 15

echo "=== DONE ==="
