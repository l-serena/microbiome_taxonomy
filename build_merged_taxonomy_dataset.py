#!/usr/bin/env python3
"""
Join DNABERT-2 embeddings with Kraken/MetaPhlAn reconcile labels.
Writes per-rank taxonomic labels for hierarchical + GNN training.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np
import pandas as pd

from taxon_utils import (
    STANDARD_RANKS,
    is_labeled,
    labels_by_rank,
    load_ktaxonomy,
    norm_read_id,
)


def load_reconcile_tsv(path: Path) -> dict:
    """read_id -> row dict"""
    out = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rid = norm_read_id(row["read_id"])
            out[rid] = row
    return out


def load_embeddings(sample: str, emb_dir: Path) -> tuple[np.ndarray, pd.DataFrame]:
    npy = emb_dir / "embeddings.npy"
    meta = emb_dir / "read_metadata.tsv"
    if not npy.exists() or not meta.exists():
        raise FileNotFoundError(f"Missing embeddings for {sample}: {emb_dir}")
    emb = np.load(npy)
    meta_df = pd.read_csv(meta, sep="\t")
    if len(emb) != len(meta_df):
        raise ValueError(f"Embedding/meta length mismatch for {sample}")
    return emb, meta_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids_file", default=os.path.expanduser("~/microbiome_taxonomy/ids.txt"))
    ap.add_argument("--output_root", default=os.path.expanduser("~/scratch/output"))
    ap.add_argument("--emb_root", default=os.path.expanduser("~/scratch/dnabert2_embeddings"))
    ap.add_argument("--ktaxonomy", default=os.path.expanduser("~/scratch/kraken_db/ktaxonomy.tsv"))
    ap.add_argument("--outdir", default=os.path.expanduser("~/scratch/merged_taxonomy_ml"))
    ap.add_argument("--max_reads_per_sample", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    tax = load_ktaxonomy(args.ktaxonomy)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    rows = []
    emb_blocks = []

    with open(args.ids_file) as fh:
        samples = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]

    for sample in samples:
        final_tsv = Path(args.output_root) / sample / "final_per_read.tsv"
        emb_dir = Path(args.emb_root) / sample

        if not final_tsv.exists():
            print(f"skip {sample}: missing {final_tsv}", flush=True)
            continue
        if not emb_dir.exists():
            print(f"skip {sample}: missing embeddings {emb_dir}", flush=True)
            continue

        labels = load_reconcile_tsv(final_tsv)
        emb, meta = load_embeddings(sample, emb_dir)

        matched_idx = []
        for i, rid in enumerate(meta["read_id"].astype(str)):
            ridn = norm_read_id(rid)
            if ridn in labels:
                matched_idx.append(i)

        if not matched_idx:
            print(f"skip {sample}: no read_id overlap", flush=True)
            continue

        if args.max_reads_per_sample and len(matched_idx) > args.max_reads_per_sample:
            matched_idx = rng.choice(matched_idx, size=args.max_reads_per_sample, replace=False).tolist()

        base = len(rows)
        emb_blocks.append(emb[matched_idx])

        for local_i, i in enumerate(matched_idx):
            rid = norm_read_id(str(meta.iloc[i]["read_id"]))
            rec = labels[rid]
            final_taxid = str(rec.get("final_taxid", "0")).strip() or "0"
            kraken_taxid = str(rec.get("kraken_filtered_taxid", rec.get("kraken_taxid", "0"))).strip() or "0"
            meta_taxid = str(rec.get("meta_taxid", "0")).strip() or "0"
            rank_labels = labels_by_rank(final_taxid, tax)

            row = {
                "sample": sample,
                "read_id": rid,
                "embedding_index": base + local_i,
                "final_taxid": final_taxid,
                "final_rank": rec.get("final_rank", ""),
                "final_name": rec.get("final_name", ""),
                "kraken_taxid": kraken_taxid,
                "meta_taxid": meta_taxid,
                "meta_hit": rec.get("meta_hit", ""),
                "reason": rec.get("reason", ""),
                "is_labeled": int(is_labeled(final_taxid)),
            }
            for r in STANDARD_RANKS:
                row[f"label_{r}"] = rank_labels[r]
            rows.append(row)

        print(f"{sample}: matched={len(matched_idx)}", flush=True)

    if not rows:
        raise SystemExit("No merged rows produced. Run pipeline + DNABERT embeddings first.")

    df = pd.DataFrame(rows)
    embeddings = np.vstack(emb_blocks).astype(np.float32)

    manifest = outdir / "merged_manifest.tsv"
    emb_path = outdir / "embeddings.npy"
    df.to_csv(manifest, sep="\t", index=False)
    np.save(emb_path, embeddings)

    # train/test split by sample (avoid leakage)
    samples_present = sorted(df["sample"].unique())
    rng_split = np.random.default_rng(args.seed)
    rng_split.shuffle(samples_present)
    n_test = max(1, int(0.2 * len(samples_present)))
    test_samples = set(samples_present[:n_test])
    train_mask = ~df["sample"].isin(test_samples)

    train_df = df[train_mask].copy().reset_index(drop=True)
    test_df = df[~train_mask].copy().reset_index(drop=True)
    train_emb = embeddings[train_mask.to_numpy()]
    test_emb = embeddings[(~train_mask).to_numpy()]

    train_df.to_csv(outdir / "train_manifest.tsv", sep="\t", index=False)
    test_df.to_csv(outdir / "test_manifest.tsv", sep="\t", index=False)
    np.save(outdir / "train_embeddings.npy", train_emb)
    np.save(outdir / "test_embeddings.npy", test_emb)

    print(f"wrote {manifest} rows={len(df)} emb={embeddings.shape}", flush=True)
    print(f"train={len(train_df)} test={len(test_df)}", flush=True)


if __name__ == "__main__":
    main()
