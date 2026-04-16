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

def load_sample_labels(final_tsv):
    labels = {}
    with open(final_tsv) as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r, None)
        for row in r:
            if len(row) < 9:
                continue
            read_id = norm_read_id(row[0])
            final_taxid = str(row[6]).strip()
            final_rank = row[7].strip()
            final_name = row[8].strip()

            if not final_taxid or final_taxid == "0":
                continue

            labels[read_id] = (final_taxid, final_rank, final_name)
    return labels

def stream_fastq_to_csv(sample, fastq_path, labels, writer):
    matched = 0
    with gzip.open(fastq_path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()

            rid = norm_read_id(h)
            if rid in labels:
                taxid, rank, name = labels[rid]
                writer.writerow([sample, rid, seq, taxid, rank, name])
                matched += 1
    return matched

def main():
    sample_dirs = sorted([d for d in OUTPUTS.iterdir() if d.is_dir()])[:10]
    if not sample_dirs:
        raise ValueError("No sample directories found in ~/scratch/outputs")

    combined_csv = str(OUT_PREFIX) + ".csv"
    train_csv = str(OUT_PREFIX) + ".train.csv"
    test_csv = str(OUT_PREFIX) + ".test.csv"

    total_rows = 0
    with open(combined_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["sample", "read_id", "sequence", "label_taxid", "label_rank", "label_name"])

        for sample_dir in sample_dirs:
            sample = sample_dir.name
            final_tsv = sample_dir / "final_per_read.tsv"

            if not final_tsv.exists():
                print("skip missing final_per_read:", sample, flush=True)
                continue

            labels = load_sample_labels(final_tsv)
            if not labels:
                print("skip no labels:", sample, flush=True)
                continue

            fq1 = DEHOST / "{}_dehost_1.fastq.gz".format(sample)
            fq2 = DEHOST / "{}_dehost_2.fastq.gz".format(sample)
            fqse = DEHOST / "{}_dehost.fastq.gz".format(sample)

            matched = 0
            if fq1.exists():
                matched += stream_fastq_to_csv(sample, fq1, labels, w)
            if fq2.exists():
                matched += stream_fastq_to_csv(sample, fq2, labels, w)
            if fqse.exists():
                matched += stream_fastq_to_csv(sample, fqse, labels, w)

            total_rows += matched
            print("{} labels={} matched_sequences={}".format(sample, len(labels), matched), flush=True)

    print("wrote combined:", combined_csv, "rows=", total_rows, flush=True)

    df = pd.read_csv(combined_csv)

    counts = df["label_taxid"].value_counts()
    keep = set(counts[counts >= 2].index.astype(str))
    df["label_taxid"] = df["label_taxid"].astype(str)
    df = df[df["label_taxid"].isin(keep)].copy()

    if len(df) == 0:
        raise ValueError("No rows left after requiring >=2 examples per class")

    train_df, test_df = train_test_split(
        df,
        test_size=0.2,
        random_state=42,
        stratify=df["label_taxid"],
    )

    train_df.to_csv(train_csv, index=False)
    test_df.to_csv(test_csv, index=False)

    print("wrote train:", train_csv, len(train_df), flush=True)
    print("wrote test:", test_csv, len(test_df), flush=True)

if __name__ == "__main__":
    main()
