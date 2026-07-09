#!/usr/bin/env python3
"""
Greedy hierarchical multiclass classifier:
train one classifier per taxonomic rank, conditioned on parent rank label.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.preprocessing import LabelEncoder

from taxon_utils import STANDARD_RANKS


def parent_rank(rank: str) -> str | None:
    if rank not in STANDARD_RANKS:
        return None
    i = STANDARD_RANKS.index(rank)
    return STANDARD_RANKS[i - 1] if i > 0 else None


def build_parent_features(df: pd.DataFrame, rank: str, encoders: dict) -> np.ndarray:
    p = parent_rank(rank)
    if p is None or p not in encoders:
        return np.zeros((len(df), 0), dtype=np.float32)
    col = f"label_{p}"
    le: LabelEncoder = encoders[p]
    parents = df[col].astype(str).fillna("0")
  # unseen -> 0
    known = set(le.classes_)
    parents = parents.where(parents.isin(known), "0")
    return le.transform(parents).reshape(-1, 1).astype(np.float32)


def train_level(
    rank: str,
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    train_emb: np.ndarray,
    test_emb: np.ndarray,
    encoders: dict,
    min_class_count: int,
):
    col = f"label_{rank}"
    train_y = train_df[col].astype(str)
    test_y = test_df[col].astype(str)

    # only labeled reads
    train_mask = train_df["is_labeled"].astype(bool).to_numpy()
    test_mask = test_df["is_labeled"].astype(bool).to_numpy()

    counts = train_y[train_mask].value_counts()
    keep = set(counts[counts >= min_class_count].index.astype(str))
    keep.add("0")

    tr_idx = train_mask & train_y.isin(keep)
    te_idx = test_mask & test_y.isin(keep)

    if tr_idx.sum() < 50:
        return None, {"rank": rank, "skipped": True, "reason": "too_few_rows"}

    le = LabelEncoder()
    le.fit(sorted(keep))
    encoders[rank] = le

    parent_tr = build_parent_features(train_df.loc[tr_idx], rank, encoders)
    parent_te = build_parent_features(test_df.loc[te_idx], rank, encoders)
    X_tr = np.hstack([train_emb[tr_idx.to_numpy()], parent_tr]) if parent_tr.size else train_emb[tr_idx.to_numpy()]
    X_te = np.hstack([test_emb[te_idx.to_numpy()], parent_te]) if parent_te.size else test_emb[te_idx.to_numpy()]
    y_tr = le.transform(train_y[tr_idx])
    y_te = le.transform(test_y[te_idx])

    clf = LogisticRegression(max_iter=2000, class_weight="balanced", n_jobs=4)
    clf.fit(X_tr, y_tr)
    pred = clf.predict(X_te)

    metrics = {
        "rank": rank,
        "n_train": int(tr_idx.sum()),
        "n_test": int(te_idx.sum()),
        "n_classes": int(len(le.classes_)),
        "accuracy": float(accuracy_score(y_te, pred)),
        "macro_f1": float(f1_score(y_te, pred, average="macro", zero_division=0)),
    }
    return clf, metrics


def greedy_predict(
    emb: np.ndarray,
    df: pd.DataFrame,
    models: dict,
    encoders: dict,
) -> pd.DataFrame:
    out = df.copy()
    for rank in STANDARD_RANKS:
        col = f"pred_{rank}"
        if rank not in models:
            out[col] = "0"
            continue
        le = encoders[rank]
        p = parent_rank(rank)
        if p is None:
            parent_feat = np.zeros((len(df), 0), dtype=np.float32)
        else:
            parent_col = f"pred_{p}"
            parents = out[parent_col].astype(str).fillna("0")
            known = set(le.classes_)
            parents = parents.where(parents.isin(known), "0")
            parent_feat = le.transform(parents).reshape(-1, 1).astype(np.float32)
        X = np.hstack([emb, parent_feat]) if parent_feat.size else emb
        pred_ids = models[rank].predict(X)
        out[col] = le.inverse_transform(pred_ids)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default=os.path.expanduser("~/scratch/merged_taxonomy_ml"))
    ap.add_argument("--outdir", default=os.path.expanduser("~/scratch/models/hierarchical_greedy"))
    ap.add_argument("--min_class_count", type=int, default=20)
    args = ap.parse_args()

    datadir = Path(args.datadir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    train_df = pd.read_csv(datadir / "train_manifest.tsv", sep="\t")
    test_df = pd.read_csv(datadir / "test_manifest.tsv", sep="\t")
    train_emb = np.load(datadir / "train_embeddings.npy")
    test_emb = np.load(datadir / "test_embeddings.npy")

    encoders: dict = {}
    models: dict = {}
    metrics = []

    for rank in STANDARD_RANKS:
        clf, m = train_level(
            rank, train_df, test_df, train_emb, test_emb, encoders, args.min_class_count
        )
        metrics.append(m)
        print(m, flush=True)
        if clf is not None:
            models[rank] = clf

    with open(outdir / "models.pkl", "wb") as fh:
        pickle.dump({"models": models, "encoders": encoders, "ranks": STANDARD_RANKS}, fh)

    pd.DataFrame(metrics).to_csv(outdir / "metrics.tsv", sep="\t", index=False)

    pred_df = greedy_predict(test_emb, test_df, models, encoders)
    pred_df.to_csv(outdir / "test_predictions.tsv", sep="\t", index=False)

    # species-level accuracy when labeled
    if "label_species" in pred_df.columns and "pred_species" in pred_df.columns:
        mask = pred_df["is_labeled"].astype(bool)
        acc = (pred_df.loc[mask, "label_species"] == pred_df.loc[mask, "pred_species"]).mean()
        print(f"species accuracy (labeled): {acc:.4f}", flush=True)


if __name__ == "__main__":
    main()
