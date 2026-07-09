#!/usr/bin/env python3
"""Shared taxonomy helpers for reconcile + ML dataset builders."""

from __future__ import annotations

STANDARD_RANKS = [
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def norm_read_id(x: str) -> str:
    x = x.strip()
    if x.startswith("@") or x.startswith(">"):
        x = x[1:]
    x = x.split()[0]
    x = x.split("__")[0]
    if x.endswith("/1") or x.endswith("/2"):
        x = x[:-2]
    return x


def load_ktaxonomy(path: str) -> dict:
    tax = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                taxid, parent, rank, name = parts[0], parts[1], parts[2], parts[3]
                tax[taxid] = {"parent": parent, "rank": rank, "name": name}
    return tax


def lineage(taxid: str, tax: dict) -> list[str]:
    out = []
    seen = set()
    t = str(taxid) if taxid else "0"
    while t and t != "0" and t not in seen:
        seen.add(t)
        out.append(t)
        t = tax.get(t, {}).get("parent", "")
    return list(reversed(out))


def rank_of(taxid: str, tax: dict) -> str:
    if not taxid or taxid == "0":
        return "unclassified"
    return tax.get(str(taxid), {}).get("rank", "unclassified")


def name_of(taxid: str, tax: dict) -> str:
    if not taxid or taxid == "0":
        return "unclassified"
    return tax.get(str(taxid), {}).get("name", "unclassified")


def labels_by_rank(final_taxid: str, tax: dict, ranks: list[str] | None = None) -> dict[str, str]:
    ranks = ranks or STANDARD_RANKS
    out = {r: "0" for r in ranks}
    for tid in lineage(final_taxid, tax):
        r = rank_of(tid, tax)
        if r in out:
            out[r] = tid
    return out


def is_labeled(final_taxid: str) -> bool:
    return bool(final_taxid) and str(final_taxid) != "0"
