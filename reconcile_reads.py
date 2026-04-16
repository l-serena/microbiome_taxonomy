#!/usr/bin/env python3
import argparse, bz2, csv, os, ast

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

# kraken taxonomy
tax = {}
with open(args.ktaxonomy) as fh:
    for l in fh:
        p = l.rstrip("\n").split("\t")
        if len(p) >= 4:
            taxid, parent, rank, name = p[0], p[1], p[2], p[3]
            tax[taxid] = {"parent": parent, "rank": rank, "name": name}

# metaphlan profile: full clade_name -> last numeric taxid in NCBI_tax_id column
clade2taxid = {}
with open(args.metaphlan_profile) as fh:
    for l in fh:
        if l.startswith("#"):
            continue
        p = l.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0] != "clade_name":
            clade2taxid[p[0]] = p[1].split("|")[-1]

# marker_info: marker -> full lineage string from dict["taxon"]
marker2clade = {}
with open_file(args.marker_info) as fh:
    for l in fh:
        p = l.rstrip("\n").split("\t", 1)
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
    for l in fh:
        if l.startswith("@"):
            continue
        p = l.rstrip("\n").split("\t")
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
        "kraken_taxid", "kraken_rank", "kraken_name",
        "meta_hit", "meta_taxid",
        "final_taxid", "final_rank", "final_name", "reason"
    ])

    with open(args.kraken_perread) as fh:
        for l in fh:
            p = l.rstrip("\n").split("\t")
            if len(p) < 3:
                continue

            status, rid, ktaxid = p[0], p[1], p[2]
            if status == "U":
                ktaxid = "0"

            kr_rank = tax.get(ktaxid, {}).get("rank", "unclassified")
            kr_name = tax.get(ktaxid, {}).get("name", "unclassified")

            ridn = norm(rid)
            meta_taxid = ""

            if ridn in read2marker:
                marker = read2marker[ridn]
                clade = marker2clade.get(marker, "")
                meta_taxid = clade2taxid.get(clade, "")

            if not meta_taxid:
                final = ktaxid
                reason = "no_meta"
            else:
                if meta_taxid in lineage(ktaxid, tax) or ktaxid in lineage(meta_taxid, tax):
                    if len(lineage(ktaxid, tax)) >= len(lineage(meta_taxid, tax)):
                        final = ktaxid
                    else:
                        final = meta_taxid
                    reason = "same_lineage"
                else:
                    final = lca(ktaxid, meta_taxid, tax)
                    reason = "lca"

            fr = tax.get(final, {}).get("rank", "unclassified")
            fn = tax.get(final, {}).get("name", "unclassified")

            w.writerow([
                rid,
                ktaxid, kr_rank, kr_name,
                "yes" if meta_taxid else "no", meta_taxid,
                final, fr, fn, reason
            ])
