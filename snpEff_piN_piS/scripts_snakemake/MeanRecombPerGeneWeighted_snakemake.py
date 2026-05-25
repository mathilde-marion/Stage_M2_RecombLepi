#!/usr/bin/env python3

"""
Compute length-weighted mean recombination statistics per gene
from bedtools intersect -wa -wb output.

Input format (9 columns):
gene_chr gene_start gene_end gene_name
win_chr  win_start  win_end  recomb_rate  recomb_quantile_100

Outputs:
- recombRate_mean: length-weighted mean recombination rate per gene
- recombQuantile_mean: length-weighted mean recombination quantile per gene
"""

import sys
from collections import defaultdict

intersect_file = sys.argv[1]
out_tsv = sys.argv[2]

# gene -> [sum(rate * overlap), sum(overlap), sum(q * overlap)]
acc = defaultdict(lambda: [0.0, 0, 0.0])

with open(intersect_file) as f:
    for line in f:
        parts = line.strip().split("\t")

        gene_start = int(parts[1])
        gene_end   = int(parts[2])
        gene       = parts[3]

        win_start  = int(parts[5])
        win_end    = int(parts[6])

        rate       = float(parts[7])
        quantile   = float(parts[8])

        overlap = max(0, min(gene_end, win_end) - max(gene_start, win_start))

        if overlap > 0:
            acc[gene][0] += rate * overlap
            acc[gene][1] += overlap
            acc[gene][2] += quantile * overlap

with open(out_tsv, "w") as out:
    out.write("gene_id\trecombRate\trecombQuantile\n")
    for gene, (sum_rate, denom, sum_q) in acc.items():
        if denom > 0:
            out.write(f"{gene}\t{sum_rate/denom}\t{sum_q/denom}\n")
