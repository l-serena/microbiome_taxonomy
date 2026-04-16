#!/bin/bash
set -euo pipefail

source /apps/lib-osver/miniforge3/etc/profile.d/conda.sh
conda activate dnabert

echo "=== Python ==="
which python
python --version

echo "=== STEP 0: package check ==="
python - <<'PY'
mods = ["pandas", "sklearn", "datasets", "transformers", "evaluate", "torch"]
missing = []
for m in mods:
    try:
        __import__(m)
    except Exception:
        missing.append(m)
if missing:
    raise SystemExit("Missing packages in dnabert env: " + ", ".join(missing))
print("All required packages found.")
PY

echo "=== STEP 1: build and split dataset ==="
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
