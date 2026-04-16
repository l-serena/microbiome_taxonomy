#!/usr/bin/env bash
set -euo pipefail

IDS_FILE="${1:?Usage: bash run_all_samples.sh ids.txt [threads]}"
THREADS="${2:-8}"

REPO="$HOME/microbiome_taxonomy"
SCRATCH="$HOME/scratch"

FASTQ_DIR="$SCRATCH/fastq"
SRA_DIR="$SCRATCH/sra"
TMP_DIR="$SCRATCH/tmp"
DEHOST_DIR="$SCRATCH/dehost_fastq"
OUT_ROOT="$SCRATCH/outputs"

HUMAN_BT2_PREFIX="$SCRATCH/human_db/GRCh38_noalt_as/GRCh38_noalt_as"
KRAKEN_DB="$SCRATCH/kraken_db"
KTAX="$KRAKEN_DB/ktaxonomy.tsv"

METAPHLAN_DB="$SCRATCH/metaphlan_db"
MPA_VER="mpa_vJan25_CHOCOPhlAnSGB_202503"
MARKER_INFO="$METAPHLAN_DB/${MPA_VER}_marker_info.txt.bz2"

RECON_SCRIPT="$REPO/reconcile_reads.py"

mkdir -p "$FASTQ_DIR" "$SRA_DIR" "$TMP_DIR" "$DEHOST_DIR" "$OUT_ROOT"

# ---------------------------
# conda init
# ---------------------------
CONDA_EXE="/apps/lib-osver/miniforge3/bin/conda"
if [[ ! -x "$CONDA_EXE" ]]; then
  echo "ERROR: conda not found at $CONDA_EXE" >&2
  exit 1
fi

__conda_setup="$("$CONDA_EXE" shell.bash hook 2>/dev/null)" || true
if [[ -n "${__conda_setup:-}" ]]; then
  eval "$__conda_setup"
else
  source /apps/lib-osver/miniforge3/etc/profile.d/conda.sh
fi
unset __conda_setup

# ---------------------------
# sanity checks
# ---------------------------
if ! ls "${HUMAN_BT2_PREFIX}".*.bt2* >/dev/null 2>&1; then
  echo "ERROR: missing Bowtie2 human index: $HUMAN_BT2_PREFIX" >&2
  exit 1
fi

if [[ ! -f "$KTAX" ]]; then
  echo "ERROR: missing ktaxonomy.tsv: $KTAX" >&2
  exit 1
fi

if [[ ! -f "$MARKER_INFO" ]]; then
  echo "ERROR: missing MetaPhlAn marker_info: $MARKER_INFO" >&2
  exit 1
fi

if [[ ! -f "$RECON_SCRIPT" ]]; then
  echo "ERROR: missing reconciliation script: $RECON_SCRIPT" >&2
  exit 1
fi

run_one() {
  local ACC="$1"

  [[ -z "$ACC" ]] && return 0
  [[ "$ACC" =~ ^[[:space:]]*$ ]] && return 0
  [[ "$ACC" =~ ^# ]] && return 0

  local SAMPLE_OUT="$OUT_ROOT/$ACC"
  mkdir -p "$SAMPLE_OUT"

  local FASTQ1="$FASTQ_DIR/${ACC}_1.fastq.gz"
  local FASTQ2="$FASTQ_DIR/${ACC}_2.fastq.gz"
  local FASTQSE="$FASTQ_DIR/${ACC}.fastq.gz"

  local DEHOST1="$DEHOST_DIR/${ACC}_dehost_1.fastq.gz"
  local DEHOST2="$DEHOST_DIR/${ACC}_dehost_2.fastq.gz"
  local DEHOSTSE="$DEHOST_DIR/${ACC}_dehost.fastq.gz"

  local KRAKEN_PERREAD="$SAMPLE_OUT/kraken2.perread.tsv"
  local KRAKEN_REPORT="$SAMPLE_OUT/kraken2.report.txt"

  local META_PROFILE="$SAMPLE_OUT/metaphlan.profile.tsv"
  local META_BOWTIE="$SAMPLE_OUT/metaphlan.bowtie2out.bz2"

  local FINAL_OUT="$SAMPLE_OUT/final_per_read.tsv"

  echo "=================================================="
  echo "Sample: $ACC"
  echo "Out:    $SAMPLE_OUT"
  echo "=================================================="

  # --------------------------------------------------
  # 1) download FASTQ if needed
  # --------------------------------------------------
  if [[ ( -f "$FASTQ1" && -f "$FASTQ2" ) || -f "$FASTQSE" ]]; then
    echo "[$ACC] FASTQ exists"
  else
    echo "[$ACC] Downloading from SRA..."
    conda activate sra

    prefetch "$ACC" --output-directory "$SRA_DIR" >"$SAMPLE_OUT/prefetch.log" 2>&1

    local SRA_PATH=""
    if [[ -f "$SRA_DIR/$ACC/$ACC.sra" ]]; then
      SRA_PATH="$SRA_DIR/$ACC/$ACC.sra"
    elif [[ -f "$SRA_DIR/$ACC.sra" ]]; then
      SRA_PATH="$SRA_DIR/$ACC.sra"
    else
      SRA_PATH="$(find "$SRA_DIR" -maxdepth 3 -name "${ACC}.sra" -print -quit || true)"
    fi

    if [[ -z "$SRA_PATH" ]]; then
      echo "ERROR: [$ACC] .sra not found" >&2
      conda deactivate || true
      return 1
    fi

    fasterq-dump "$SRA_PATH" \
      --outdir "$FASTQ_DIR" \
      --temp "$TMP_DIR" \
      --threads "$THREADS" \
      --split-files \
      >"$SAMPLE_OUT/fasterq.log" 2>&1

    conda deactivate || true

    if [[ -f "$FASTQ_DIR/${ACC}_1.fastq" && -f "$FASTQ_DIR/${ACC}_2.fastq" ]]; then
      gzip -f "$FASTQ_DIR/${ACC}_1.fastq" "$FASTQ_DIR/${ACC}_2.fastq"
    elif [[ -f "$FASTQ_DIR/${ACC}.fastq" ]]; then
      gzip -f "$FASTQ_DIR/${ACC}.fastq"
    else
      echo "ERROR: [$ACC] no FASTQ produced" >&2
      return 1
    fi
  fi

  # --------------------------------------------------
  # 2) dehost
  # --------------------------------------------------
  if [[ ( -f "$DEHOST1" && -f "$DEHOST2" ) || -f "$DEHOSTSE" ]]; then
    echo "[$ACC] dehosted FASTQ exists"
  else
    echo "[$ACC] Host removal..."
    conda activate biobakery3

    if [[ -f "$FASTQ1" && -f "$FASTQ2" ]]; then
      bowtie2 \
        -x "$HUMAN_BT2_PREFIX" \
        -1 "$FASTQ1" -2 "$FASTQ2" \
        -p "$THREADS" \
        --very-sensitive \
        --un-conc-gz "$DEHOST_DIR/${ACC}_dehost_%.fastq.gz" \
        -S /dev/null \
        2>"$SAMPLE_OUT/bowtie2.log"
    else
      local INSE=""
      [[ -f "$FASTQSE" ]] && INSE="$FASTQSE"
      [[ -z "$INSE" && -f "$FASTQ1" ]] && INSE="$FASTQ1"

      if [[ -z "$INSE" ]]; then
        echo "ERROR: [$ACC] no FASTQ input for dehosting" >&2
        conda deactivate || true
        return 1
      fi

      bowtie2 \
        -x "$HUMAN_BT2_PREFIX" \
        -U "$INSE" \
        -p "$THREADS" \
        --very-sensitive \
        --un-gz "$DEHOSTSE" \
        -S /dev/null \
        2>"$SAMPLE_OUT/bowtie2.log"
    fi

    conda deactivate || true
  fi

  # --------------------------------------------------
  # 3) kraken2
  # --------------------------------------------------
  if [[ -f "$KRAKEN_PERREAD" && -f "$KRAKEN_REPORT" ]]; then
    echo "[$ACC] Kraken2 output exists"
  else
    echo "[$ACC] Running Kraken2..."
    conda activate kraken2

    if [[ -f "$DEHOST1" && -f "$DEHOST2" ]]; then
      kraken2 \
        --db "$KRAKEN_DB" \
        --threads "$THREADS" \
        --paired \
        --confidence 0.1 \
        --report "$KRAKEN_REPORT" \
        --output "$KRAKEN_PERREAD" \
        "$DEHOST1" "$DEHOST2" \
        >"$SAMPLE_OUT/kraken2.log" 2>&1
    else
      local INSE=""
      [[ -f "$DEHOSTSE" ]] && INSE="$DEHOSTSE"
      [[ -z "$INSE" && -f "$DEHOST1" ]] && INSE="$DEHOST1"

      if [[ -z "$INSE" ]]; then
        echo "ERROR: [$ACC] no dehosted FASTQ for Kraken2" >&2
        conda deactivate || true
        return 1
      fi

      kraken2 \
        --db "$KRAKEN_DB" \
        --threads "$THREADS" \
        --confidence 0.1 \
        --report "$KRAKEN_REPORT" \
        --output "$KRAKEN_PERREAD" \
        "$INSE" \
        >"$SAMPLE_OUT/kraken2.log" 2>&1
    fi

    conda deactivate || true
  fi

  # --------------------------------------------------
  # 4) MetaPhlAn
  # --------------------------------------------------
  if [[ -f "$META_PROFILE" && -f "$META_BOWTIE" ]]; then
    echo "[$ACC] MetaPhlAn output exists"
  else
    echo "[$ACC] Running MetaPhlAn..."
    conda activate biobakery3

    if [[ -f "$DEHOST1" && -f "$DEHOST2" ]]; then
      metaphlan \
        "$DEHOST1","$DEHOST2" \
        --input_type fastq \
        --nproc "$THREADS" \
        --bowtie2db "$METAPHLAN_DB" \
        --bowtie2out "$META_BOWTIE" \
        -o "$META_PROFILE" \
        >"$SAMPLE_OUT/metaphlan.log" 2>&1
    else
      local INSE=""
      [[ -f "$DEHOSTSE" ]] && INSE="$DEHOSTSE"
      [[ -z "$INSE" && -f "$DEHOST1" ]] && INSE="$DEHOST1"

      if [[ -z "$INSE" ]]; then
        echo "ERROR: [$ACC] no dehosted FASTQ for MetaPhlAn" >&2
        conda deactivate || true
        return 1
      fi

      metaphlan \
        "$INSE" \
        --input_type fastq \
        --nproc "$THREADS" \
        --bowtie2db "$METAPHLAN_DB" \
        --bowtie2out "$META_BOWTIE" \
        -o "$META_PROFILE" \
        >"$SAMPLE_OUT/metaphlan.log" 2>&1
    fi

    conda deactivate || true
  fi

  # --------------------------------------------------
  # 5) reconcile
  # --------------------------------------------------
  if [[ -f "$FINAL_OUT" ]]; then
    echo "[$ACC] final_per_read.tsv exists"
  else
    echo "[$ACC] Reconciling..."
    conda activate biobakery3

    python "$RECON_SCRIPT" \
      --kraken_perread "$KRAKEN_PERREAD" \
      --metaphlan_bowtie2 "$META_BOWTIE" \
      --metaphlan_profile "$META_PROFILE" \
      --marker_info "$MARKER_INFO" \
      --ktaxonomy "$KTAX" \
      --out "$FINAL_OUT" \
      >"$SAMPLE_OUT/reconcile.log" 2>&1

    conda deactivate || true
  fi

  echo "[$ACC] done"
}

while IFS= read -r ACC || [[ -n "$ACC" ]]; do
  ACC="$(echo "$ACC" | tr -d '\r')"
  run_one "$ACC"
done < "$IDS_FILE"

echo "All samples complete."
