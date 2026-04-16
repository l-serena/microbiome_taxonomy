#!/bin/bash
#SBATCH --job-name=tax_pipeline
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail

IDS_FILE=$1

ACC=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$IDS_FILE" | tr -d '\r')

if [[ -z "$ACC" || "$ACC" =~ ^[[:space:]]*$ || "$ACC" =~ ^# ]]; then
  echo "Skipping empty line"
  exit 0
fi

THREADS=${SLURM_CPUS_PER_TASK:-8}

BASE=$HOME/scratch
OUT=$BASE/outputs/$ACC
mkdir -p "$OUT"

DEHOST=$BASE/dehost_fastq
KDB=$BASE/kraken_db
MDB=$BASE/metaphlan_db
KTAX=$KDB/ktaxonomy.tsv
MARKER=$MDB/mpa_vJan25_CHOCOPhlAnSGB_202503_marker_info.txt.bz2
REPO=$HOME/microbiome_taxonomy

KRAKEN_PERREAD="$OUT/kraken2.perread.tsv"
KRAKEN_REPORT="$OUT/kraken2.report.txt"
META_BOWTIE="$OUT/metaphlan.bowtie2out.bz2"
META_PROFILE="$OUT/metaphlan.profile.tsv"
FINAL_OUT="$OUT/final_per_read.tsv"

R1="$DEHOST/${ACC}_dehost_1.fastq.gz"
R2="$DEHOST/${ACC}_dehost_2.fastq.gz"
SE="$DEHOST/${ACC}_dehost.fastq.gz"

echo "Running $ACC"

source /apps/lib-osver/miniforge3/etc/profile.d/conda.sh

conda activate kraken2

if [[ ! -s "$KRAKEN_PERREAD" || ! -s "$KRAKEN_REPORT" ]]; then
  echo "[$ACC] Running Kraken2..."

  if [[ -f "$R1" && -f "$R2" ]]; then
    kraken2 \
      --db "$KDB" \
      --threads "$THREADS" \
      --paired \
      --confidence 0.1 \
      --report "$KRAKEN_REPORT" \
      --output "$KRAKEN_PERREAD" \
      "$R1" "$R2" \
      > "$OUT/kraken.log" 2>&1
  elif [[ -f "$SE" ]]; then
    kraken2 \
      --db "$KDB" \
      --threads "$THREADS" \
      --confidence 0.1 \
      --report "$KRAKEN_REPORT" \
      --output "$KRAKEN_PERREAD" \
      "$SE" \
      > "$OUT/kraken.log" 2>&1
  else
    echo "ERROR: no dehosted FASTQ found for $ACC" >&2
    exit 1
  fi
fi

if [[ ! -s "$KRAKEN_PERREAD" || ! -s "$KRAKEN_REPORT" ]]; then
  echo "ERROR: Kraken2 failed for $ACC" >&2
  echo "See $OUT/kraken.log" >&2
  exit 1
fi

conda deactivate
conda activate biobakery3

if [[ ! -s "$META_PROFILE" || ! -s "$META_BOWTIE" ]]; then
  echo "[$ACC] Running MetaPhlAn..."

  if [[ -f "$R1" && -f "$R2" ]]; then
    metaphlan \
      "$R1","$R2" \
      --input_type fastq \
      --nproc "$THREADS" \
      --bowtie2db "$MDB" \
      --bowtie2out "$META_BOWTIE" \
      -o "$META_PROFILE" \
      > "$OUT/metaphlan.log" 2>&1
  elif [[ -f "$SE" ]]; then
    metaphlan \
      "$SE" \
      --input_type fastq \
      --nproc "$THREADS" \
      --bowtie2db "$MDB" \
      --bowtie2out "$META_BOWTIE" \
      -o "$META_PROFILE" \
      > "$OUT/metaphlan.log" 2>&1
  else
    echo "ERROR: no dehosted FASTQ found for MetaPhlAn: $ACC" >&2
    exit 1
  fi
fi

if [[ ! -s "$META_PROFILE" || ! -s "$META_BOWTIE" ]]; then
  echo "ERROR: MetaPhlAn failed for $ACC" >&2
  echo "See $OUT/metaphlan.log" >&2
  exit 1
fi

rm -f "$FINAL_OUT"

python "$REPO/reconcile_reads.py" \
  --kraken_perread "$KRAKEN_PERREAD" \
  --metaphlan_bowtie2 "$META_BOWTIE" \
  --metaphlan_profile "$META_PROFILE" \
  --marker_info "$MARKER" \
  --ktaxonomy "$KTAX" \
  --out "$FINAL_OUT" \
  > "$OUT/reconcile.log" 2>&1

if [[ ! -s "$FINAL_OUT" ]]; then
  echo "ERROR: reconciliation failed for $ACC" >&2
  echo "See $OUT/reconcile.log" >&2
  exit 1
fi

conda deactivate

echo "DONE $ACC"
