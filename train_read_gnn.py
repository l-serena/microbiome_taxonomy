#!/usr/bin/env python3
"""
Read-level GNN classifier using DNABERT embeddings + kNN graph within each sample.
Unknown/unlabeled reads are masked in the loss.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import f1_score
from sklearn.neighbors import NearestNeighbors

from taxon_utils import STANDARD_RANKS


class ReadGNN(nn.Module):
    def __init__(self, in_dim: int, hidden: int, n_classes: int):
        super().__init__()
        self.lin1 = nn.Linear(in_dim, hidden)
        self.lin2 = nn.Linear(hidden, hidden)
        self.out = nn.Linear(hidden, n_classes)

    def forward(self, x, adj):
        # adj: row-normalized sparse or dense [N,N]
        h = F.relu(self.lin1(x))
        h = torch.matmul(adj, h)
        h = F.relu(self.lin2(h))
        h = torch.matmul(adj, h)
        return self.out(h)


def build_knn_adj(x: np.ndarray, k: int) -> np.ndarray:
    n = x.shape[0]
    k = min(k, max(1, n - 1))
    nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine")
    nn.fit(x)
    dist, idx = nn.kneighbors(x)
    adj = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in idx[i]:
            if i != j:
                adj[i, j] = 1.0
        adj[i, i] = 1.0
    deg = adj.sum(axis=1, keepdims=True)
    adj = adj / np.clip(deg, 1e-6, None)
    return adj


def train_epoch(model, optimizer, samples, device):
    model.train()
    total_loss = 0.0
    n_batches = 0
    for emb, y, mask in samples:
        emb_t = torch.tensor(emb, dtype=torch.float32, device=device)
        y_t = torch.tensor(y, dtype=torch.long, device=device)
        mask_t = torch.tensor(mask, dtype=torch.bool, device=device)
        adj = torch.tensor(build_knn_adj(emb, k=8), dtype=torch.float32, device=device)

        logits = model(emb_t, adj)
        if mask_t.sum() == 0:
            continue
        loss = F.cross_entropy(logits[mask_t], y_t[mask_t])
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        total_loss += float(loss.item())
        n_batches += 1
    return total_loss / max(n_batches, 1)


@torch.no_grad()
def eval_samples(model, samples, device):
    model.eval()
    preds, refs, masks = [], [], []
    for emb, y, mask in samples:
        emb_t = torch.tensor(emb, dtype=torch.float32, device=device)
        adj = torch.tensor(build_knn_adj(emb, k=8), dtype=torch.float32, device=device)
        logits = model(emb_t, adj)
        pred = logits.argmax(dim=-1).cpu().numpy()
        preds.append(pred)
        refs.append(y)
        masks.append(mask)
    preds = np.concatenate(preds)
    refs = np.concatenate(refs)
    masks = np.concatenate(masks).astype(bool)
    if masks.sum() == 0:
        return {"accuracy": 0.0, "macro_f1": 0.0}
    acc = (preds[masks] == refs[masks]).mean()
    f1 = f1_score(refs[masks], preds[masks], average="macro", zero_division=0)
    return {"accuracy": float(acc), "macro_f1": float(f1)}


def group_by_sample(df: pd.DataFrame, emb: np.ndarray, label_col: str, le_classes: dict):
    assert len(df) == len(emb)
    groups = []
    for sample, sub in df.groupby("sample"):
        idx = sub.index.to_numpy()
        x = emb[idx]
        y_raw = sub[label_col].astype(str)
        y = np.array([le_classes.get(v, -1) for v in y_raw], dtype=np.int64)
        mask = (sub["is_labeled"].astype(bool).to_numpy()) & (y >= 0)
        groups.append((x, y, mask))
    return groups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default=os.path.expanduser("~/scratch/merged_taxonomy_ml"))
    ap.add_argument("--outdir", default=os.path.expanduser("~/scratch/models/read_gnn"))
    ap.add_argument("--label_rank", default="species", choices=STANDARD_RANKS)
    ap.add_argument("--min_class_count", type=int, default=20)
    ap.add_argument("--hidden", type=int, default=256)
    ap.add_argument("--epochs", type=int, default=15)
    ap.add_argument("--lr", type=float, default=1e-3)
    args = ap.parse_args()

    datadir = Path(args.datadir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    train_df = pd.read_csv(datadir / "train_manifest.tsv", sep="\t")
    test_df = pd.read_csv(datadir / "test_manifest.tsv", sep="\t")
    train_emb = np.load(datadir / "train_embeddings.npy")
    test_emb = np.load(datadir / "test_embeddings.npy")

    label_col = f"label_{args.label_rank}"
    counts = train_df.loc[train_df["is_labeled"].astype(bool), label_col].astype(str).value_counts()
    keep = sorted(counts[counts >= args.min_class_count].index.astype(str))
    le_classes = {c: i for i, c in enumerate(keep)}

    tr_mask = train_df[label_col].astype(str).isin(keep)
    train_sub = train_df[tr_mask].reset_index(drop=True)
    train_emb_sub = train_emb[tr_mask.to_numpy()]

    te_mask = test_df[label_col].astype(str).isin(keep)
    test_sub = test_df[te_mask].reset_index(drop=True)
    test_emb_sub = test_emb[te_mask.to_numpy()]

    train_groups = group_by_sample(train_sub, train_emb_sub, label_col, le_classes)
    test_groups = group_by_sample(test_sub, test_emb_sub, label_col, le_classes)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = ReadGNN(train_emb.shape[1], args.hidden, len(keep)).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    best_f1 = -1.0
    for epoch in range(1, args.epochs + 1):
        loss = train_epoch(model, optimizer, train_groups, device)
        metrics = eval_samples(model, test_groups, device)
        print(f"epoch {epoch} loss={loss:.4f} {metrics}", flush=True)
        if metrics["macro_f1"] > best_f1:
            best_f1 = metrics["macro_f1"]
            torch.save(
                {
                    "state_dict": model.state_dict(),
                    "label_rank": args.label_rank,
                    "classes": keep,
                    "hidden": args.hidden,
                },
                outdir / "best_model.pt",
            )

    pd.DataFrame([{"label_rank": args.label_rank, "best_macro_f1": best_f1}]).to_csv(
        outdir / "metrics.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
