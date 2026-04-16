#!/usr/bin/env python3
import argparse, ast, bz2, csv, os

def open_file(p):
    return bz2.open(p, "rt") if p.endswith(".bz2") else open(p)

def norm(x):
    x = x.strip()
    if x.startswith("@") or x.startswith(">"):
        x = x[1:]
    x = x.split()[0]
    x = x.split("__")[0]
    if x.endswith("/1") or x.endswith("/2"):
        x = x[:-2]
    return x

def lineage(t, tax):
    out = []
    seen = set()
    while t and t != "0" and t not in seen:
        seen.add(t)
        out.append(t)
        t = tax.get(t, {}).get("parent", "")
    return list(reversed(out))

def lca(a, b, tax):
    A, B = lineage(a, tax), lineage(b, tax)
    x = "0"
    for i, j in zip(A, B):
        if i == j:
            x = i
        else:
            break
    return x

def rank_of(t, tax):
    if not t or t == "0":
        return "unclassified"
    return tax.get(t, {}).get("rank", "unclassified")

def name_of(t, tax):
    if not t or t == "0":
        return "unclassified"
    return tax.get(t, {}).get("name", "unclassified")

def ancestor_at_rank(t, wanted_rank, tax):
    while t and t != "0":
        if rank_of(t, tax) == wanted_rank:
            return t
        t = tax.get(t, {}).get("parent", "")
    return "0"

def kraken_confidence(kmer_field, target_taxid):
    total = 0
    support = 0
    target_lineage = set()

    if target_taxid and target_taxid != "0":
        target_lineage = set(lineage(target_taxid, TAX))

    for token in kmer_field.split():
        if token == "|:|":
            continue
        if ":" not in token:
            continue
        tid, count = token.split(":", 1)
        try:
            count = int(count)
        except Exception:
            continue
        total += count
        if tid in target_lineage:
            support += count

    return (support / total) if total > 0 else 0.0

ap = argparse.ArgumentParser()
ap.add_argument("--kraken_perread", required=True)
ap.add_argument("--metaphlan_bowtie2", required=True)
ap.add_argument("--metaphlan_profile", required=True)
ap.add_argument("--marker_info", required=True)
ap.add_argument("--ktaxonomy", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

for p in [args.kraken_perread, args.metaphlan_bowtie2, args.metaphlan_profile, args.marker_info, args.ktaxonomy]:
    if not os.path.exists(p):
        raise FileNotFoundError(p)

# Kraken taxonomy
TAX = {}
with open(args.ktaxonomy) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            taxid, parent, rank, name = p[0], p[1], p[2], p[3]
            TAX[taxid] = {"parent": parent, "rank": rank, "name": name}

# MetaPhlAn profile: full clade -> last non-empty taxid
clade2taxid = {}
with open(args.metaphlan_profile) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0] != "clade_name":
            ids = [x for x in p[1].split("|") if x]
            clade2taxid[p[0]] = ids[-1] if ids else ""

# marker_info: marker -> full lineage string
marker2clade = {}
with open_file(args.marker_info) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t", 1)
        if len(p) < 2:
            continue
        marker, info = p[0], p[1]
        try:
            d = ast.literal_eval(info)
            clade = d.get("taxon", "")
            if clade:
                marker2clade[marker] = clade
        except Exception:
            pass

# bowtie2out: read -> marker
read2marker = {}
with open_file(args.metaphlan_bowtie2) as fh:
    for line in fh:
        if line.startswith("@"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) >= 11:
            rid, marker = p[0], p[2]
        elif len(p) >= 2:
            rid, marker = p[0], p[1]
        else:
            continue
        rid = norm(rid)
        if rid not in read2marker:
            read2marker[rid] = marker

with open(args.out, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow([
        "read_id",
        "kraken_raw_taxid", "kraken_raw_rank", "kraken_raw_name", "kraken_confidence",
        "kraken_filtered_taxid", "kraken_filtered_rank", "kraken_filtered_name",
        "meta_hit", "meta_taxid", "meta_rank", "meta_name",
        "final_taxid", "final_rank", "final_name", "reason"
    ])

    with open(args.kraken_perread) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue

            status, rid, raw_taxid = p[0], p[1], p[2]
            kmer_field = p[4] if len(p) >= 5 else ""

            if status == "U":
                raw_taxid = "0"

            raw_rank = rank_of(raw_taxid, TAX)
            raw_name = name_of(raw_taxid, TAX)
            conf = kraken_confidence(kmer_field, raw_taxid)

            # Kraken backoff logic
            if raw_taxid == "0":
                k_taxid = "0"
            elif raw_rank == "species" and conf >= 0.2:
                k_taxid = raw_taxid
            else:
                genus_taxid = ancestor_at_rank(raw_taxid, "genus", TAX)
                genus_conf = kraken_confidence(kmer_field, genus_taxid) if genus_taxid != "0" else 0.0

                if genus_taxid != "0" and genus_conf >= 0.1:
                    k_taxid = genus_taxid
                else:
                    family_taxid = ancestor_at_rank(raw_taxid, "family", TAX)
                    k_taxid = family_taxid if family_taxid != "0" else raw_taxid

            k_rank = rank_of(k_taxid, TAX)
            k_name = name_of(k_taxid, TAX)

            ridn = norm(rid)
            meta_taxid = ""
            if ridn in read2marker:
                marker = read2marker[ridn]
                clade = marker2clade.get(marker, "")
                meta_taxid = clade2taxid.get(clade, "")

            meta_rank = rank_of(meta_taxid, TAX) if meta_taxid else ""
            meta_name = name_of(meta_taxid, TAX) if meta_taxid else ""

            # Reconciliation
            if not meta_taxid:
                final = k_taxid
                reason = "kraken_only"
            else:
                if k_taxid == "0":
                    final = meta_taxid
                    reason = "meta_only"
                elif k_taxid == meta_taxid:
                    final = k_taxid
                    reason = "exact_match"
                elif meta_taxid in lineage(k_taxid, TAX) or k_taxid in lineage(meta_taxid, TAX):
                    final = k_taxid if len(lineage(k_taxid, TAX)) >= len(lineage(meta_taxid, TAX)) else meta_taxid
                    reason = "same_lineage"
                else:
                    final = lca(k_taxid, meta_taxid, TAX)
                    reason = "lca"

            final_rank = rank_of(final, TAX)
            final_name = name_of(final, TAX)

            w.writerow([
                rid,
                raw_taxid, raw_rank, raw_name, f"{conf:.6f}",
                k_taxid, k_rank, k_name,
                "yes" if meta_taxid else "no", meta_taxid, meta_rank, meta_name,
                final, final_rank, final_name, reason
            ])
