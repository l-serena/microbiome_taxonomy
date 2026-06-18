#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# run_one_sample.sh (SLURM array worker)
# ============================================================

IDS_FILE="${1:?Usage: $0 ids.txt}"

ACC="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$IDS_FILE" | tr -d '\r' | tr -d '[:space:]')"
if [[ -z "${ACC}" || "${ACC}" =~ ^[[:space:]]*$ || "${ACC}" =~ ^# ]]; then
  echo "Task ${SLURM_ARRAY_TASK_ID}: empty/comment line; exiting."
  exit 0
fi

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ---- SLURM-safe conda init ----
CONDA_EXE="/apps/lib-osver/miniforge3/bin/conda"
if [[ ! -x "$CONDA_EXE" ]]; then
  echo "ERROR: conda executable not found at $CONDA_EXE" >&2
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
# Paths
# ---------------------------
BASE="$HOME/scratch"
FASTQ_DIR="$BASE/fastq"
SRA_DIR="$BASE/sra"
TMP_DIR="${SLURM_TMPDIR:-$BASE/tmp}"
DEHOST_DIR="$BASE/dehost_fastq"

OUTPUT_ROOT="$BASE/output"
SAMPLE_OUT="$OUTPUT_ROOT/$ACC"

mkdir -p "$FASTQ_DIR" "$SRA_DIR" "$TMP_DIR" "$DEHOST_DIR" "$SAMPLE_OUT/logs"

# ---------------------------
# DBs / Indexes
# ---------------------------
HUMAN_BT2_INDEX_PREFIX="$BASE/human_db/GRCh38_noalt_as/GRCh38_noalt_as"
KRAKEN_DB="$BASE/kraken_db"
METAPHLAN_DB="$BASE/metaphlan_db"

# Marker info for your MetaPhlAn DB
MPA_VER="mpa_vJan25_CHOCOPhlAnSGB_202503"
MARKER_INFO="$METAPHLAN_DB/${MPA_VER}_marker_info.txt.bz2"
MARKER_INFO_URL="http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/${MPA_VER}_marker_info.txt.bz2"

# ---------------------------
# Inputs/Outputs
# ---------------------------
FASTQ1="$FASTQ_DIR/${ACC}_1.fastq.gz"
FASTQ2="$FASTQ_DIR/${ACC}_2.fastq.gz"
FASTQSE="$FASTQ_DIR/${ACC}.fastq.gz"

DEHOST1="$DEHOST_DIR/${ACC}_dehost_1.fastq.gz"
DEHOST2="$DEHOST_DIR/${ACC}_dehost_2.fastq.gz"
DEHOSTSE="$DEHOST_DIR/${ACC}_dehost.fastq.gz"

KRAKEN_PERREAD="$SAMPLE_OUT/kraken2.perread.tsv"
KRAKEN_REPORT="$SAMPLE_OUT/kraken2.report.txt"

META_PROFILE="$SAMPLE_OUT/metaphlan.profile.tsv"
META_BOWTIE="$SAMPLE_OUT/metaphlan.bowtie2out.bz2"

FINAL_DIR="$SAMPLE_OUT/reconcile"
FINAL_OUTDIR="$FINAL_DIR/out"
FINAL_OUT_SENTINEL="$FINAL_OUTDIR/final_per_read"   # directory created by reconciler

echo "===================================="
echo "ACC=$ACC"
echo "THREADS=$THREADS"
echo "HOSTNAME=$(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-na}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-na}"
echo "FASTQ_DIR=$FASTQ_DIR"
echo "DEHOST_DIR=$DEHOST_DIR"
echo "SAMPLE_OUT=$SAMPLE_OUT"
echo "===================================="

# ------------------------------------------------------------
# Sanity: Bowtie2 human index exists
# ------------------------------------------------------------
if ! ls "${HUMAN_BT2_INDEX_PREFIX}".*.bt2* >/dev/null 2>&1; then
  echo "ERROR: Bowtie2 index not found for prefix: $HUMAN_BT2_INDEX_PREFIX" >&2
  exit 1
fi

# ------------------------------------------------------------
# 1) FASTQ CHECK
#    FASTQs should already exist in $HOME/scratch/fastq
#    Accepts:
#      paired-end: <ACC>_1.fastq.gz and <ACC>_2.fastq.gz
#      single-end: <ACC>.fastq.gz
# ------------------------------------------------------------
if [[ -f "$FASTQ1" && -f "$FASTQ2" ]]; then
  echo "==> [$ACC] Paired-end FASTQs found."
elif [[ -f "$FASTQSE" ]]; then
  echo "==> [$ACC] Single-end FASTQ found."
else
  echo "ERROR: [$ACC] FASTQ files not found in $FASTQ_DIR" >&2
  echo "Expected either:" >&2
  echo "  $FASTQ1 and $FASTQ2" >&2
  echo "or:" >&2
  echo "  $FASTQSE" >&2
  exit 1
fi

# ------------------------------------------------------------
# 2) HOST REMOVAL
# ------------------------------------------------------------
if [[ ( -f "$DEHOST1" && -f "$DEHOST2" ) || -f "$DEHOSTSE" ]]; then
  echo "==> [$ACC] Dehosted reads exist — skipping bowtie2."
else
  echo "==> [$ACC] Host removal with bowtie2..."
  conda activate biobakery3

  if [[ -f "$FASTQ1" && -f "$FASTQ2" ]]; then
    bowtie2 -x "$HUMAN_BT2_INDEX_PREFIX" \
      -1 "$FASTQ1" -2 "$FASTQ2" \
      -p "$THREADS" --very-sensitive \
      --un-conc-gz "$DEHOST_DIR/${ACC}_dehost_%.fastq.gz" \
      -S /dev/null 2>"$SAMPLE_OUT/logs/${ACC}.bowtie2.log"
  else
    INSE=""
    [[ -f "$FASTQSE" ]] && INSE="$FASTQSE"
    [[ -z "$INSE" && -f "$FASTQ1" ]] && INSE="$FASTQ1"
    if [[ -z "$INSE" ]]; then
      echo "ERROR: [$ACC] No FASTQ input found for host removal." >&2
      conda deactivate || true
      exit 1
    fi
    bowtie2 -x "$HUMAN_BT2_INDEX_PREFIX" \
      -U "$INSE" \
      -p "$THREADS" --very-sensitive \
      --un-gz "$DEHOSTSE" \
      -S /dev/null 2>"$SAMPLE_OUT/logs/${ACC}.bowtie2.log"
  fi

  conda deactivate || true
fi

# ------------------------------------------------------------
# 3) KRAKEN2
# ------------------------------------------------------------
if [[ -f "$KRAKEN_PERREAD" && -f "$KRAKEN_REPORT" ]]; then
  echo "==> [$ACC] Kraken outputs exist — skipping."
else
  echo "==> [$ACC] Running Kraken2..."
  conda activate kraken2

  if [[ -f "$DEHOST1" && -f "$DEHOST2" ]]; then
    kraken2 --db "$KRAKEN_DB" --threads "$THREADS" --paired --confidence 0.1 \
      --report "$KRAKEN_REPORT" --output "$KRAKEN_PERREAD" \
      "$DEHOST1" "$DEHOST2" \
      >"$SAMPLE_OUT/logs/${ACC}.kraken2.log" 2>&1
  else
    INSE=""
    [[ -f "$DEHOSTSE" ]] && INSE="$DEHOSTSE"
    [[ -z "$INSE" && -f "$DEHOST1" ]] && INSE="$DEHOST1"
    if [[ -z "$INSE" ]]; then
      echo "ERROR: [$ACC] No dehost FASTQ input found for Kraken2." >&2
      conda deactivate || true
      exit 1
    fi
    kraken2 --db "$KRAKEN_DB" --threads "$THREADS" --confidence 0.1 \
      --report "$KRAKEN_REPORT" --output "$KRAKEN_PERREAD" \
      "$INSE" \
      >"$SAMPLE_OUT/logs/${ACC}.kraken2.log" 2>&1
  fi

  conda deactivate || true
fi

# ------------------------------------------------------------
# 4) METAPHLAN (resumable)
# ------------------------------------------------------------
if [[ -f "$META_PROFILE" ]]; then
  echo "==> [$ACC] MetaPhlAn profile exists — skipping."
elif [[ -f "$META_BOWTIE" ]]; then
  echo "==> [$ACC] Bowtie2out exists — generating MetaPhlAn profile from bowtie2out..."
  conda activate biobakery3
  metaphlan "$META_BOWTIE" --input_type bowtie2out --nproc "$THREADS" --bowtie2db "$METAPHLAN_DB" -o "$META_PROFILE" \
    >"$SAMPLE_OUT/logs/${ACC}.metaphlan.profile_only.log" 2>&1
  conda deactivate || true
else
  echo "==> [$ACC] Running MetaPhlAn (mapping + profile)..."
  conda activate biobakery3
  if [[ -f "$DEHOST1" && -f "$DEHOST2" ]]; then
    metaphlan "$DEHOST1","$DEHOST2" --input_type fastq --nproc "$THREADS" --bowtie2db "$METAPHLAN_DB" --bowtie2out "$META_BOWTIE" -o "$META_PROFILE" \
      >"$SAMPLE_OUT/logs/${ACC}.metaphlan.log" 2>&1
  else
    INSE=""
    [[ -f "$DEHOSTSE" ]] && INSE="$DEHOSTSE"
    [[ -z "$INSE" && -f "$DEHOST1" ]] && INSE="$DEHOST1"
    if [[ -z "$INSE" ]]; then
      echo "ERROR: [$ACC] No dehost FASTQ input found for MetaPhlAn." >&2
      conda deactivate || true
      exit 1
    fi
    metaphlan "$INSE" --input_type fastq --nproc "$THREADS" --bowtie2db "$METAPHLAN_DB" --bowtie2out "$META_BOWTIE" -o "$META_PROFILE" \
      >"$SAMPLE_OUT/logs/${ACC}.metaphlan.log" 2>&1
  fi
  conda deactivate || true
fi

# ------------------------------------------------------------
# 4.5) Ensure marker_info exists (auto-download via HTTP)
# ------------------------------------------------------------
if [[ -f "$MARKER_INFO" ]]; then
  echo "==> [$ACC] marker_info exists — ok."
else
  echo "==> [$ACC] marker_info missing; downloading (HTTP)..."
  mkdir -p "$METAPHLAN_DB"
  wget -O "$MARKER_INFO" "$MARKER_INFO_URL" \
    >"$SAMPLE_OUT/logs/${ACC}.marker_info.download.log" 2>&1 || {
      echo "ERROR: failed to download marker_info from $MARKER_INFO_URL" >&2
      echo "See $SAMPLE_OUT/logs/${ACC}.marker_info.download.log" >&2
      exit 1
    }
fi

# ------------------------------------------------------------
# 5) RECONCILE (directory-based, Python 3 in biobakery3)
#    Expects:
#      metaphlan_dir: <ACC>_1.bowtie2.bz2 and <ACC>_1.profile.tsv
#      kraken_dir:    <ACC>.kraken
# ------------------------------------------------------------
if [[ -d "$FINAL_OUT_SENTINEL" ]]; then
  echo "==> [$ACC] Reconciliation outputs exist — skipping."
else
  echo "==> [$ACC] Reconciling per-read labels..."
  mkdir -p "$FINAL_DIR/metaphlan_dir" "$FINAL_DIR/kraken_dir" "$FINAL_OUTDIR"

  ln -sf "$META_BOWTIE"  "$FINAL_DIR/metaphlan_dir/${ACC}_1.bowtie2.bz2"
  ln -sf "$META_PROFILE" "$FINAL_DIR/metaphlan_dir/${ACC}_1.profile.tsv"
  ln -sf "$KRAKEN_PERREAD" "$FINAL_DIR/kraken_dir/${ACC}.kraken"

  conda activate biobakery3
  python -u reconcile_reads_special.py \
    --metaphlan_dir "$FINAL_DIR/metaphlan_dir" \
    --kraken_dir "$FINAL_DIR/kraken_dir" \
    --kraken_db_root "$KRAKEN_DB" \
    --marker_info "$MARKER_INFO" \
    --outdir "$FINAL_OUTDIR" \
    >"$SAMPLE_OUT/logs/${ACC}.reconcile.log" 2>&1
  conda deactivate || true
fi

echo "==> [$ACC] DONE"
echo "Sample outputs:"
echo "  Kraken:   $KRAKEN_REPORT"
echo "  KrakenPR: $KRAKEN_PERREAD"
echo "  Meta:     $META_PROFILE"
echo "  MetaB2:   $META_BOWTIE"
echo "  Reconcile:$FINAL_OUTDIR"
echo "Logs:       $SAMPLE_OUT/logs/"
