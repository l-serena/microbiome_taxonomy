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
        header = next(r, None)
        for row in r:
            if len(row) < 9:
                continue

            read_id = norm_read_id(row[0])
            final_taxid = str(row[6]).strip()
            final_rank = row[7].strip()
            final_name = row[8].strip()

            if not final_taxid or final_taxid == "0":
                continue

            labels[read_id] = {
                "taxid": final_taxid,
                "rank": final_rank,
                "name": final_name,
            }
    return labels

def main():
    sample_dirs = sorted([d for d in OUTPUTS.iterdir() if d.is_dir()])[:10]
    if not sample_dirs:
        raise ValueError("No sample directories found in ~/scratch/outputs")

    rows = []

    for sample_dir in sample_dirs:
        sample = sample_dir.name
        final_tsv = sample_dir / "final_per_read.tsv"

        if not final_tsv.exists():
            print("skip missing final_per_read:", sample)
            continue

        labels = load_sample_labels(final_tsv)
        if not labels:
            print("skip no usable labels:", sample)
            continue

        fq1 = DEHOST / "{}_dehost_1.fastq.gz".format(sample)
        fq2 = DEHOST / "{}_dehost_2.fastq.gz".format(sample)
        fqse = DEHOST / "{}_dehost.fastq.gz".format(sample)

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

        print("{} labels={} matched_sequences={}".format(sample, len(labels), matched))

    df = pd.DataFrame(rows)
    if len(df) == 0:
        raise ValueError("No matched read/label pairs found")

    # Need at least 2 per class for stratified split
    counts = df["label_taxid"].value_counts()
    keep = set(counts[counts >= 2].index)
    df = df[df["label_taxid"].isin(keep)].copy()

    if len(df) == 0:
        raise ValueError("No rows left after requiring >=2 examples per class")

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

    print("wrote {} rows to {}".format(len(train_df), train_path))
    print("wrote {} rows to {}".format(len(test_df), test_path))

if __name__ == "__main__":
    main()
