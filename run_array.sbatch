#!/bin/bash
#SBATCH --job-name=tax_pipeline
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00

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

echo "Running $ACC"

source /apps/lib-osver/miniforge3/etc/profile.d/conda.sh

conda activate kraken2
if [[ ! -f "$OUT/kraken2.perread.tsv" ]]; then
  kraken2 \
    --db "$KDB" \
    --threads "$THREADS" \
    --paired \
    --confidence 0.1 \
    --report "$OUT/kraken2.report.txt" \
    --output "$OUT/kraken2.perread.tsv" \
    "$DEHOST/${ACC}_dehost_1.fastq.gz" \
    "$DEHOST/${ACC}_dehost_2.fastq.gz" \
    > "$OUT/kraken.log" 2>&1
fi
conda deactivate

conda activate biobakery3
if [[ ! -f "$OUT/metaphlan.profile.tsv" ]]; then
  metaphlan \
    "$DEHOST/${ACC}_dehost_1.fastq.gz","$DEHOST/${ACC}_dehost_2.fastq.gz" \
    --input_type fastq \
    --nproc "$THREADS" \
    --bowtie2db "$MDB" \
    --bowtie2out "$OUT/metaphlan.bowtie2out.bz2" \
    -o "$OUT/metaphlan.profile.tsv" \
    > "$OUT/metaphlan.log" 2>&1
fi

if [[ ! -f "$OUT/final_per_read.tsv" ]]; then
  python "$REPO/reconcile_reads.py" \
    --kraken_perread "$OUT/kraken2.perread.tsv" \
    --metaphlan_bowtie2 "$OUT/metaphlan.bowtie2out.bz2" \
    --metaphlan_profile "$OUT/metaphlan.profile.tsv" \
    --marker_info "$MARKER" \
    --ktaxonomy "$KTAX" \
    --out "$OUT/final_per_read.tsv" \
    > "$OUT/reconcile.log" 2>&1
fi
conda deactivate

echo "DONE $ACC"
