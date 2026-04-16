#!/usr/bin/env python3
import csv
import gzip
import os
from pathlib import Path

OUTPUTS = Path(os.path.expanduser("~/scratch/outputs"))
DEHOST = Path(os.path.expanduser("~/scratch/dehost_fastq"))
OUT_CSV = Path(os.path.expanduser("~/scratch/dnabert2_reads_10samples.csv"))

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
            rid = norm_read_id(h)
            seqs[rid] = s
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
                "taxid": final_taxid,
                "rank": final_rank,
                "name": final_name,
            }
    return labels

def main():
    sample_dirs = sorted([d for d in OUTPUTS.iterdir() if d.is_dir()])[:10]

    rows_written = 0
    with open(OUT_CSV, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["sample", "read_id", "sequence", "label_taxid", "label_rank", "label_name"])

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
                    w.writerow([sample, rid, seqs[rid], lab["taxid"], lab["rank"], lab["name"]])
                    matched += 1
                    rows_written += 1

            print(sample, "labels=", len(labels), "matched_sequences=", matched)

    print("wrote:", OUT_CSV)
    print("total rows:", rows_written)

if __name__ == "__main__":
    main()
