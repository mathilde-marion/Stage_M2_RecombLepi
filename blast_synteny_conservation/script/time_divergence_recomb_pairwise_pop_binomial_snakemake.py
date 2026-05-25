#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute proportion of conserved recombination rates (success)
per species pair and test significance using a binomial test,
then add divergence time from a dated phylogenetic tree.

Input:
  - TSV with columns: sp1, sp2, success (0/1)
  - Dated phylogenetic tree (Newick)

Output:
  - TSV with:
      sp1, sp2,
      n_windows, n_success, mean_success,
      pvalue_binom,
      divergence_time_Ma
"""

import argparse
import pandas as pd
from Bio import Phylo
from scipy.stats import binomtest

# =========================
# ARGUMENTS
# =========================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Add divergence time to pairwise recombination statistics"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input TSV file with recombination statistics"
    )

    parser.add_argument(
        "-t", "--tree",
        required=True,
        help="Dated phylogenetic tree in Newick format"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file"
    )

    return parser.parse_args()

# =========================
# SPECIES NAME MAPPING
# =========================

NAME_MAP = {
    "ito_sal": "Itomia_salapia",
    "lys_bel": "Lysandra_bellargus",
    "pol_ica": "Polyommatus_icarus",
    "mel_pho": "Melitaea_phoebe",
    "pie_nap": "Pieris_napi",
    "pap_mac": "Papilio_machaon",
    "pap_hos": "Papilio_hospiton",
    "hel_mel": "Heliconius_melpomene",
    "hel_tim": "Heliconius_timareta",
    "hel_cyd": "Heliconius_cydno",
    "hel_num": "Heliconius_numata",
    "ant_car": "Anthocharis_cardamines",
    "mel_gal": "Melanargia_galathea"
}

def sp_code(sp_pop):
    """Extrait le code espèce depuis 'hel_mel_POP1' → 'hel_mel'"""
    # Suppose que le suffixe pop est toujours _POP1 ou _POP2
    return "_".join(sp_pop.split("_")[:-1])

# =========================
# MAIN
# =========================

def main():

    args = parse_args()

    # ---- Load recombination data
    df = pd.read_csv(args.input, sep="\t")

    if "success" not in df.columns:
        raise ValueError("❌ Column 'success' not found in input file")

    # ---- Compute binomial stats per species pair
    def compute_stats(group):
        n = len(group)
        k = group["success"].sum()
        mean_success = k / n if n > 0 else 0

        res = binomtest(k, n, p=0.5, alternative="greater")

        return pd.Series({
            "n_windows": n,
            "n_success": k,
            "mean_success": mean_success,
            "pvalue_binom": res.pvalue
        })

    pairwise = (
        df
        .groupby(["sp1", "sp2"], as_index=False)
        .apply(compute_stats)
        .reset_index(drop=True)
    )

    # ---- Map short names to tree names
#    pairwise["sp1_tree"] = pairwise["sp1"].map(NAME_MAP)
#    pairwise["sp2_tree"] = pairwise["sp2"].map(NAME_MAP)
#
#    if pairwise[["sp1_tree", "sp2_tree"]].isna().any().any():
#        raise ValueError("Some species could not be mapped to tree names")

    # ---- Map short names to tree names
    pairwise["sp1_tree"] = pairwise["sp1"].apply(lambda x: NAME_MAP.get(sp_code(x)))
    pairwise["sp2_tree"] = pairwise["sp2"].apply(lambda x: NAME_MAP.get(sp_code(x)))

    if pairwise[["sp1_tree", "sp2_tree"]].isna().any().any():
        missing = pairwise[pairwise[["sp1_tree","sp2_tree"]].isna().any(axis=1)][["sp1","sp2"]]
        raise ValueError(f"Species non mappées :\n{missing}")

    # ---- Load tree
    tree = Phylo.read(args.tree, "newick")
    tree_tips = {t.name for t in tree.get_terminals()}

    for sp in pd.concat([pairwise["sp1_tree"], pairwise["sp2_tree"]]):
        if sp not in tree_tips:
            raise ValueError(f"Species '{sp}' not found in tree")

    # ---- Divergence time function
    def divergence_time(tree, sp1_tree, sp2_tree):
        if sp1_tree == sp2_tree:
            return 0.0 # même espèce, pops différentes donc divergence de zéro
        clade1 = tree.find_any(name=sp1_tree)
        clade2 = tree.find_any(name=sp2_tree)
        mrca = tree.common_ancestor(clade1, clade2)
        return tree.distance(mrca, clade1)

    # ---- Compute divergence times
    pairwise["divergence_time_Ma"] = pairwise.apply(
        lambda r: divergence_time(tree, r["sp1_tree"], r["sp2_tree"]),
        axis=1
    )

    # ---- Final output
    final_df = pairwise[[
        "sp1",
        "sp2",
        "n_windows",
        "n_success",
        "mean_success",
        "pvalue_binom",
        "divergence_time_Ma"
    ]]

    final_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
