#!/usr/bin/env python3
import glob
import sys
import os
import re

blast_dir = sys.argv[1]
output_dir = sys.argv[2]
species = sys.argv[3]

os.makedirs(output_dir, exist_ok=True)

sp_dir = f"{blast_dir}/{species}"
out_path = os.path.join(output_dir, f"{species}.blast.best.tsv")

results = []

for f in glob.glob(f"{sp_dir}/*.blast.tsv"):

    gene_full = os.path.basename(f).replace(".blast.tsv", "")
    gene = gene_full.split(".")[0]   # enlève .ENSCAHG...

    best_line = None
    best_evalue = float("inf")

    with open(f) as in_f:
        for line in in_f:
            if not line.strip():
                continue

            parts = line.strip().split()
            evalue = float(parts[10])

            if evalue < best_evalue:
                best_evalue = evalue
                best_line = parts

    if best_line is None:
        best_line = [gene] + ["NA"] * 11
    else:
        best_line[0] = gene

    results.append(best_line)


# fonction de tri naturel (chr_gene)
def sort_key(row):
    g = row[0]
    m = re.match(r"(\d+)_(\d+)", g)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    return (999,999)


results.sort(key=sort_key)


with open(out_path, "w") as out_f:

    header = [
        "qseqid","sseqid","pident","length","mismatch","gapopen",
        "qstart","qend","sstart","send","evalue","bitscore"
    ]
    out_f.write("\t".join(header) + "\n")

    for row in results:
        out_f.write("\t".join(row) + "\n")


print(f"[{species}] Best hits written to {out_path}")
