#!/usr/bin/env bash
# Run on cluster after MetaPhlAn DB install completes.
# Usage:
#   bash run_full_cluster_pipeline.sh
# Or step-by-step as jobs finish.

set -euo pipefail

REPO="${REPO:-$HOME/microbiome_taxonomy}"
IDS="$REPO/ids.txt"
N=$(grep -cve '^\s*$' -e '^\s*#' "$IDS" || true)
ARRAY_MAX=$((N - 1))

echo "Samples in ids.txt: $N (array 0-$ARRAY_MAX)"

echo "=== 1) Kraken + MetaPhlAn + reconcile ==="
cd "$REPO"
sbatch --array=0-"$ARRAY_MAX" submit.sbatch "$IDS"

echo "=== 2) After step 1 completes, DNABERT embeddings ==="
echo "  sbatch --array=0-$ARRAY_MAX run_dnabert_embeddings_array.sbatch $IDS"

echo "=== 3) After step 2 completes, merged dataset + train ==="
echo "  sbatch run_ml_pipeline.sbatch"

echo "=== Monitor ==="
echo "  bash $REPO/check_pipeline_status.sh $IDS"
