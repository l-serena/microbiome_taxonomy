#!/usr/bin/env python3
import argparse
import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer


def read_fastq(path, max_reads=None):
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rt") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break

            seq = f.readline().strip()
            f.readline()
            f.readline()

            read_id = header.split()[0].lstrip("@")
            yield read_id, seq

            n += 1
            if max_reads is not None and n >= max_reads:
                break


def mean_pool(last_hidden_state, attention_mask):
    mask = attention_mask.unsqueeze(-1).float()
    summed = (last_hidden_state * mask).sum(dim=1)
    counts = mask.sum(dim=1).clamp(min=1e-9)
    return summed / counts


def get_hidden_state(model_output):
    if hasattr(model_output, "last_hidden_state"):
        return model_output.last_hidden_state

    if isinstance(model_output, (tuple, list)):
        return model_output[0]

    raise TypeError(f"Unexpected model output type: {type(model_output)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fastq_dir", default=os.path.expanduser("~/scratch/fastq"))
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--model_name", default="zhihan1996/DNABERT-2-117M")
    ap.add_argument("--max_length", type=int, default=250)
    ap.add_argument("--batch_size", type=int, default=8)
    ap.add_argument("--max_reads", type=int, default=None)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fastq_dir = Path(args.fastq_dir)
    r1 = fastq_dir / f"{args.sample}_1.fastq.gz"
    r2 = fastq_dir / f"{args.sample}_2.fastq.gz"
    se = fastq_dir / f"{args.sample}.fastq.gz"

    if r1.exists() and r2.exists():
        files = [r1, r2]
    elif se.exists():
        files = [se]
    else:
        raise FileNotFoundError(f"No FASTQ found for {args.sample} in {fastq_dir}")

    print(f"Sample: {args.sample}")
    print(f"FASTQ files: {[str(x) for x in files]}")
    print(f"Loading model: {args.model_name}")

    tokenizer = AutoTokenizer.from_pretrained(args.model_name, trust_remote_code=True)
    model = AutoModel.from_pretrained(args.model_name, trust_remote_code=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    model.to(device)
    model.eval()

    read_ids = []
    seqs = []
    mates = []

    for fq in files:
        mate = "SE"
        if fq.name.endswith("_1.fastq.gz"):
            mate = "R1"
        elif fq.name.endswith("_2.fastq.gz"):
            mate = "R2"

        for rid, seq in read_fastq(fq, args.max_reads):
            read_ids.append(rid)
            seqs.append(seq)
            mates.append(mate)

    if len(seqs) == 0:
        raise ValueError(f"No reads found for {args.sample}")

    all_emb = []

    for start in range(0, len(seqs), args.batch_size):
        batch = seqs[start:start + args.batch_size]

        toks = tokenizer(
            batch,
            truncation=True,
            padding=True,
            max_length=args.max_length,
            return_tensors="pt",
        )
        toks = {k: v.to(device) for k, v in toks.items()}

        with torch.no_grad():
            out = model(**toks, return_dict=True)
            hidden = get_hidden_state(out)
            emb = mean_pool(hidden, toks["attention_mask"])

        all_emb.append(emb.cpu().numpy())

        if start % (args.batch_size * 100) == 0:
            print(f"Processed {start}/{len(seqs)} reads", flush=True)

    embeddings = np.vstack(all_emb)

    metadata = pd.DataFrame({
        "sample": args.sample,
        "read_id": read_ids,
        "mate": mates,
        "embedding_index": np.arange(len(read_ids)),
    })

    np.save(outdir / "embeddings.npy", embeddings)
    metadata.to_csv(outdir / "read_metadata.tsv", sep="\t", index=False)

    print(f"Saved {outdir / 'embeddings.npy'}")
    print(f"Shape: {embeddings.shape}")
    print(f"Saved {outdir / 'read_metadata.tsv'}")


if __name__ == "__main__":
    main()
