#!/usr/bin/env bash
# Pipeline status for Kraken / MetaPhlAn / reconcile / DNABERT / merged ML data.
set -euo pipefail

IDS_FILE="${1:-$HOME/microbiome_taxonomy/ids.txt}"
OUT_ROOT="${OUT_ROOT:-$HOME/scratch/output}"
EMB_ROOT="${EMB_ROOT:-$HOME/scratch/dnabert2_embeddings}"
MPA_DB="${MPA_DB:-$HOME/scratch/metaphlan_db}"

echo "=== MetaPhlAn DB ==="
ls -lh "$MPA_DB"/*.pkl 2>/dev/null | head -3 || echo "no pkl"

echo
echo "=== Per-sample status ==="
printf "%-14s %-3s %-3s %-3s %-3s %-3s\n" "SAMPLE" "K" "M" "R" "E" "N"

k_ok=m_ok=r_ok=e_ok=0
n=0

while IFS= read -r ACC || [[ -n "${ACC:-}" ]]; do
  [[ -z "$ACC" || "$ACC" =~ ^# ]] && continue
  n=$((n + 1))
  S="$OUT_ROOT/$ACC"
  K="n"; M="n"; R="n"; E="n"
  [[ -s "$S/kraken2.perread.tsv" && -s "$S/kraken2.report.txt" ]] && K="y" && k_ok=$((k_ok + 1))
  [[ -s "$S/metaphlan.profile.tsv" ]] && M="y" && m_ok=$((m_ok + 1))
  [[ -s "$S/final_per_read.tsv" ]] && R="y" && r_ok=$((r_ok + 1))
  [[ -s "$EMB_ROOT/$ACC/embeddings.npy" ]] && E="y" && e_ok=$((e_ok + 1))
  printf "%-14s %-3s %-3s %-3s %-3s %-3s\n" "$ACC" "$K" "$M" "$R" "$E" "$n"
done < "$IDS_FILE"

echo
echo "Summary: samples=$n kraken=$k_ok metaphlan=$m_ok reconcile=$r_ok embeddings=$e_ok"
