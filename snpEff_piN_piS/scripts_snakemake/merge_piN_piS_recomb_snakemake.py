#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

# --------------------------------------------------
# Load inputs
# --------------------------------------------------
pi_file      = sys.argv[1]   # piN piS
sites_file   = sys.argv[2]   # syn_sites nonsyn_sites
recomb_file  = sys.argv[3]   # recombination per gene
out_file     = sys.argv[4]

pi = pd.read_csv(pi_file, sep="\t")
sites = pd.read_csv(sites_file, sep="\t")
recomb = pd.read_csv(recomb_file, sep="\t")

# --------------------------------------------------
# Fix recomb column name
# --------------------------------------------------
#recomb = recomb.rename(columns={"gene": "gene_id"})

# --------------------------------------------------
# Clean sites file (keep only valid genes)
# --------------------------------------------------
sites = sites[
    (sites["ignored_ambiguous_codons"] == 0) &
    (sites["has_terminal_stop"] == 1)
]

# --------------------------------------------------
# Merge pi + sites
# --------------------------------------------------
df = pd.merge(pi, sites, on="gene_id")

# --------------------------------------------------
# Merge recombination
# --------------------------------------------------
df = pd.merge(df, recomb, on="gene_id", how="left")

# --------------------------------------------------
# Compute normalized values safely
# --------------------------------------------------
# Use numpy to avoid division warnings and inf

df["piN_norm"] = np.divide(
    df["piN"],
    df["nonsyn_sites"],
    out=np.full(len(df), np.nan),
    where=df["nonsyn_sites"] != 0
)

df["piS_norm"] = np.divide(
    df["piS"],
    df["syn_sites"],
    out=np.full(len(df), np.nan),
    where=df["syn_sites"] != 0
)

# --------------------------------------------------
# Compute final ratio (dN/dS-like)
# --------------------------------------------------
df["ratio_piN_piS"] = np.divide(
    df["piN_norm"],
    df["piS_norm"],
    out=np.full(len(df), np.nan),
    where=df["piS_norm"] != 0
)

# --------------------------------------------------
# Output
# --------------------------------------------------
df[[
    "gene_id",
    "piN",
    "piS",
    "nonsyn_sites",
    "syn_sites",
    "piN_norm",
    "piS_norm",
    "ratio_piN_piS",
    "recombRate",
    "recombQuantile"
]].to_csv(out_file, sep="\t", index=False)
