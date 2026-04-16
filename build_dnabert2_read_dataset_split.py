#!/usr/bin/env python3
import csv
import gzip
import os
from pathlib import Path
import pandas as pd
from sklearn.model_selection import train_test_split

OUTPUTS = Path(os.path.expanduser("~/scratch/outputs"))
DEHOST = Path(os.path.expanduser("~/scratch/dehost_fastq"))
OUT_PREFIX = Path(os.path.expanduser("~/scratch/dnabert2_reads_10samples"))

def norm_read_id(x):
    x = x.strip()
    if x.startswith("@"):
        x = x[1:]
    x = x.split()[0]
    x = x.split("__")[0]
    if x.endswith("/1") or x.endswith("/2"):
        x = x[:-2]
    return x

def read_fastq_sequences(path):
    seqs = {}
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip()
            fh.readline()
            fh.readline()
            seqs[norm_read_id(h)] = s
    return seqs

def load_sample_labels(final_tsv):
    labels = {}
    with open(final_tsv) as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        for row in r:
            if len(row) < 10:
                continue
            read_id = norm_read_id(row[0])
            final_taxid = row[6]
            final_rank = row[7]
            final_name = row[8]

            if not final_taxid or final_taxid == "0":
                continue

            labels[read_id] = {
                "taxid": str(final_taxid),
                "rank": final_rank,
                "name": final_name,
            }
    return labels

def main():
    sample_dirs = sorted([d for d in OUTPUTS.iterdir() if d.is_dir()])[:10]

    rows = []
    for sample_dir in sample_dirs:
        sample = sample_dir.name
        final_tsv = sample_dir / "final_per_read.tsv"
        if not final_tsv.exists():
            print("skip missing final_per_read:", sample)
            continue

        labels = load_sample_labels(final_tsv)
        if not labels:
            print("skip no labels:", sample)
            continue

        fq1 = DEHOST / f"{sample}_dehost_1.fastq.gz"
        fq2 = DEHOST / f"{sample}_dehost_2.fastq.gz"
        fqse = DEHOST / f"{sample}_dehost.fastq.gz"

        seqs = {}
        if fq1.exists():
            seqs.update(read_fastq_sequences(fq1))
        if fq2.exists():
            seqs.update(read_fastq_sequences(fq2))
        if fqse.exists():
            seqs.update(read_fastq_sequences(fqse))

        matched = 0
        for rid, lab in labels.items():
            if rid in seqs:
                rows.append({
                    "sample": sample,
                    "read_id": rid,
                    "sequence": seqs[rid],
                    "label_taxid": lab["taxid"],
                    "label_rank": lab["rank"],
                    "label_name": lab["name"],
                })
                matched += 1

        print(sample, "labels=", len(labels), "matched_sequences=", matched)

    df = pd.DataFrame(rows)
    if len(df) == 0:
        raise ValueError("No matched rows found")

    # keep only labels with at least 2 examples so stratified split works
    counts = df["label_taxid"].value_counts()
    keep = set(counts[counts >= 2].index)
    df = df[df["label_taxid"].isin(keep)].copy()

    train_df, test_df = train_test_split(
        df,
        test_size=0.2,
        random_state=42,
        stratify=df["label_taxid"],
    )

    train_path = str(OUT_PREFIX) + ".train.csv"
    test_path = str(OUT_PREFIX) + ".test.csv"

    train_df.to_csv(train_path, index=False)
    test_df.to_csv(test_path, index=False)

    print("wrote:", train_path, len(train_df))
    print("wrote:", test_path, len(test_df))

if __name__ == "__main__":
    main()
