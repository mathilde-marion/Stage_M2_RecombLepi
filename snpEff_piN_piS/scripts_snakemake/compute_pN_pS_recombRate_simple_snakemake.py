#!/usr/bin/env python3
"""
Compute pN/pS per gene from SnpSift-filtered VCF files,
with optional recombination rate annotation.

Usage:
    python compute_pN_pS.py <syn.vcf[.gz]> <nonsyn.vcf[.gz]> <out.tsv> [recomb.tsv]

recomb.tsv format (tab-separated, with header):
    gene    recomb_rate
    GeneA   0.0023
    GeneB   0.0041
"""

import sys
import gzip
from collections import defaultdict


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def count_per_gene(vcf_path):
    """
    Count variants per gene, deduplicated by (CHROM, POS).
    A single variant annotated on N transcripts of the same gene counts as 1.
    """
    counts   = defaultdict(int)
    seen_pos = defaultdict(set)  # gene -> set of (chrom, pos) already counted

    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos   = fields[1]
            info  = fields[7]

            genes_this_variant = set()

            for entry in info.split(";"):
                if not entry.startswith("ANN="):
                    continue
                for ann in entry.replace("ANN=", "").split(","):
                    ann_fields = ann.split("|")
                    if len(ann_fields) > 3:
                        gene = ann_fields[3].strip()
                        if gene:
                            genes_this_variant.add(gene)

            for gene in genes_this_variant:
                if (chrom, pos) not in seen_pos[gene]:
                    counts[gene] += 1
                    seen_pos[gene].add((chrom, pos))

    return counts


def load_recomb(path):
    """Load a gene -> recomb_rate TSV (with header)."""
    recomb = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        next(f, None)  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                try:
                    recomb[parts[0]] = float(parts[1])
                except ValueError:
                    recomb[parts[0]] = "NA"
    return recomb


def main():
    if len(sys.argv) < 4:
        sys.exit(
            "Usage: compute_pN_pS.py <syn.vcf[.gz]> <nonsyn.vcf[.gz]> "
            "<out.tsv> [recomb.tsv]"
        )

    syn_vcf    = sys.argv[1]
    nonsyn_vcf = sys.argv[2]
    out_tsv    = sys.argv[3]
    recomb_tsv = sys.argv[4] if len(sys.argv) > 4 else None

    syn_counts    = count_per_gene(syn_vcf)
    nonsyn_counts = count_per_gene(nonsyn_vcf)

    recomb = load_recomb(recomb_tsv) if recomb_tsv else {}

    genes = sorted(set(syn_counts) | set(nonsyn_counts))

    with open(out_tsv, "w") as out:
        header = ["gene", "pS", "pN", "pN_pS"]
        if recomb:
            header.append("recombRate")
        out.write("\t".join(header) + "\n")

        for gene in genes:
            pS    = syn_counts.get(gene, 0)
            pN    = nonsyn_counts.get(gene, 0)
            ratio = f"{pN / pS:.4f}" if pS > 0 else "NA"

            row = [gene, str(pS), str(pN), ratio]
            if recomb:
                row.append(str(recomb.get(gene, "NA")))

            out.write("\t".join(row) + "\n")

    print(f"Done -> {out_tsv}  ({len(genes)} genes)")


if __name__ == "__main__":
    main()
