#!/usr/bin/env python3

"""
Compute piN and piS per gene from a SnpEff annotated VCF.
p = AC / AN (no external tool needed)

Usage:
    python compute_piN_piS.py <ann.vcf[.gz]> <summary_out.tsv>

Output:
    summary_out.tsv  : gene_id / piN / piS / piN_piS  (one row per gene)
"""

import sys
import gzip
import re
from collections import defaultdict


# open VCF (gz or plain)
def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


# extract p = AC / AN from INFO
def get_p(info):
    """
    Parse AC and AN from INFO field.
    Returns p = AC / AN, or None if missing / AN == 0.
    Multi-allelic: takes only the first AC value.
    """
    ac, an = None, None

    for field in info.split(";"):
        if field.startswith("AC="):
            ac = int(field[3:].split(",")[0])
        elif field.startswith("AN="):
            an = int(field[3:])

    if ac is None or an is None or an == 0:
        return None

    return ac / an


# parse ANN field
def parse_ann(info):
    for field in info.split(";"):
        if field.startswith("ANN="):
            return field[4:].split(",")
    return []


# classify variant as N or S
def classify_effect(ann_list):
    """
    Priority rule:
      - any missense_variant  → N
      - any synonymous_variant → S
      - otherwise             → None
    """
    is_syn = False

    for ann in ann_list:
        if "missense_variant" in ann:
            return "N"
        if "synonymous_variant" in ann:
            is_syn = True

    return "S" if is_syn else None


# extract gene_id from ANN
def get_genes(ann_list):
    """
    Return all unique gene names found in ANN (field index 3).

    Also cleans SnpEff format:
        gene:ENSG... → ENSG...
    """
    genes = set()

    for ann in ann_list:
        parts = ann.split("|")

        if len(parts) > 4 and parts[4].strip():
            gene_id = parts[4].strip()
            genes.add(gene_id)

    return genes


# natural sorting key (g1 < g2 < g10)
def natural_key(gene):
    """
    Extract numeric part for correct sorting:
        g1, g2, g10 → 1,2,10
    """
    match = re.findall(r"\d+", gene)
    return int(match[0]) if match else gene


# main
def main():

    if len(sys.argv) < 3:
        sys.exit(
            "Usage: compute_piN_piS.py <ann.vcf[.gz]> <summary_out.tsv>"
        )

    vcf_file    = sys.argv[1]
    summary_out = sys.argv[2]

    piN = defaultdict(float)
    piS = defaultdict(float)

    # deduplication: gene -> set of (chrom, pos) already seen
    seen_pos = defaultdict(set)

    with open_vcf(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            if len(cols) < 8:
                continue

            chrom = cols[0]
            pos   = cols[1]
            info  = cols[7]

            # allele frequency
            p = get_p(info)
            if p is None:
                continue

            pi = 2 * p * (1 - p)

            # annotation
            ann_list = parse_ann(info)
            effect   = classify_effect(ann_list)
            if effect is None:
                continue

            genes = get_genes(ann_list)
            if not genes:
                continue

            # accumulate per gene (deduplicated by position)
            for gene in genes:

                if (chrom, pos) in seen_pos[gene]:
                    continue
                seen_pos[gene].add((chrom, pos))

                if effect == "N":
                    piN[gene] += pi
                else:
                    piS[gene] += pi

    # write summary output
#    genes_all = sorted(set(piN) | set(piS), key=natural_key)
#    genes_all = sorted(set(map(str, piN)) | set(map(str, piS)), key=natural_key)
    genes_all = sorted(set(map(str, piN)) | set(map(str, piS)))
    
    with open(summary_out, "w") as out:
        out.write("gene_id\tpiN\tpiS\tpiN_piS\n")

        for gene in genes_all:

            pn = piN.get(gene, 0.0)
            ps = piS.get(gene, 0.0)

            ratio = pn / ps if ps > 0 else "NA"

            out.write(f"{gene}\t{pn:.6f}\t{ps:.6f}\t{ratio}\n")

    print(f"Done -> {summary_out} ({len(genes_all)} genes)")


if __name__ == "__main__":
    main()
