#!/usr/bin/env python3
import argparse, bz2, csv

def open_file(p):
    return bz2.open(p, "rt") if p.endswith(".bz2") else open(p)

def norm(x):
    x = x.strip()
    if x.startswith("@") or x.startswith(">"): x = x[1:]
    x = x.split()[0]
    x = x.split("__")[0]
    if x.endswith("/1") or x.endswith("/2"): x = x[:-2]
    return x

def lineage(t, tax):
    L = []
    while t and t != "0":
        L.append(t)
        t = tax.get(t, {}).get("parent")
    return list(reversed(L))

def lca(a,b,tax):
    A,B = lineage(a,tax), lineage(b,tax)
    out="0"
    for x,y in zip(A,B):
        if x==y: out=x
        else: break
    return out

# ---------- load data ----------

def load_tax(path):
    d={}
    for l in open(path):
        p=l.strip().split("\t")
        if len(p)>=4:
            d[p[0]]={"parent":p[1],"rank":p[2],"name":p[3]}
    return d

def load_profile(path):
    # clade -> taxid
    m={}
    for l in open(path):
        if l.startswith("#"): continue
        p=l.strip().split("\t")
        if len(p)<2 or p[0]=="clade_name": continue
        m[p[0]] = p[1].split("|")[-1]
    return m

def load_marker_info(path):
    # marker -> clade string (simple heuristic)
    m={}
    for l in open_file(path):
        p=l.strip().split("\t")
        if not p: continue
        marker=p[0]
        for x in p[1:]:
            if "k__" in x or "d__" in x:
                m[marker]=x
                break
    return m

def load_metaphlan(path):
    # read -> marker
    m={}
    for l in open_file(path):
        if l.startswith("@"): continue
        p=l.strip().split("\t")
        if len(p)>=11:
            rid,marker=p[0],p[2]
        elif len(p)>=2:
            rid,marker=p[0],p[1]
        else:
            continue
        rid=norm(rid)
        if rid not in m:
            m[rid]=marker
    return m

# ---------- main ----------

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--kraken_perread",required=True)
    ap.add_argument("--metaphlan_bowtie2",required=True)
    ap.add_argument("--metaphlan_profile",required=True)
    ap.add_argument("--marker_info",required=True)
    ap.add_argument("--ktaxonomy",required=True)
    ap.add_argument("--out",required=True)
    args=ap.parse_args()

    tax = load_tax(args.ktaxonomy)
    prof = load_profile(args.metaphlan_profile)
    marker2clade = load_marker_info(args.marker_info)
    read2marker = load_metaphlan(args.metaphlan_bowtie2)

    out=open(args.out,"w")
    w=csv.writer(out,delimiter="\t")

    w.writerow([
        "read_id",
        "kraken_taxid","kraken_rank","kraken_name",
        "meta_hit","meta_taxid",
        "final_taxid","final_rank","final_name","reason"
    ])

    for l in open(args.kraken_perread):
        p=l.strip().split("\t")
        if len(p)<3: continue

        status, rid, ktaxid = p[0], p[1], p[2]
        if status=="U": ktaxid="0"

        kr_rank = tax.get(ktaxid,{}).get("rank","unclassified")
        kr_name = tax.get(ktaxid,{}).get("name","unclassified")

        ridn = norm(rid)

        meta_taxid=""
        if ridn in read2marker:
            marker = read2marker[ridn]
            clade  = marker2clade.get(marker,"")
            meta_taxid = prof.get(clade,"")

        # ---------- reconciliation ----------
        if not meta_taxid:
            final = ktaxid
            reason = "no_meta"
        else:
            if meta_taxid in lineage(ktaxid,tax) or ktaxid in lineage(meta_taxid,tax):
                # same lineage → deeper
                final = ktaxid if len(lineage(ktaxid,tax))>=len(lineage(meta_taxid,tax)) else meta_taxid
                reason = "same_lineage"
            else:
                final = lca(ktaxid,meta_taxid,tax)
                reason = "lca"

        fr = tax.get(final,{}).get("rank","unclassified")
        fn = tax.get(final,{}).get("name","unclassified")

        w.writerow([rid, ktaxid, kr_rank, kr_name,
                    "yes" if meta_taxid else "no", meta_taxid,
                    final, fr, fn, reason])

    out.close()


if __name__=="__main__":
    main()
