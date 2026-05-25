#!/usr/bin/env python3

import pandas as pd
import argparse
import os


def key_from_pi(path):
    # pap_hos_POP1.windowed.pi -> pap_hos_POP1
    return os.path.basename(path).replace(".windowed.pi", "")


def key_from_recomb(path):
    # pie_nap.var....POP1.quantiles.tsv -> pie_nap_POP1
    base = os.path.basename(path)
    base = base.replace(".quantiles.tsv", "")
    parts = base.split(".")
    # species + POP1
    return f"{parts[0]}_{parts[-1]}"


def main(pi_files, recomb_files, output):

    pi_map = {}
    recomb_map = {}

    # --- PI ---
    for f in pi_files:
        key = key_from_pi(f)
        df = pd.read_csv(f, sep="\t")
        pi_map[key] = df["PI"].mean()

    # --- RECOMB ---
    for f in recomb_files:
        key = key_from_recomb(f)
        df = pd.read_csv(f, sep="\t")
        df.columns = [c.strip() for c in df.columns]

        rec = pd.to_numeric(df["recombRate"], errors="coerce")

        recomb_map[key] = {
            "mean": rec.mean(),
            "var": rec.var()
        }

    # --- MERGE ---
    rows = []
    for k in sorted(pi_map.keys()):
        rows.append({
            "species_pop": k,
            "pi_mean": pi_map[k],
            "recomb_mean": recomb_map[k]["mean"],
            "recomb_var": recomb_map[k]["var"]
        })

    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--pi", nargs="+", required=True)
    parser.add_argument("--recomb", nargs="+", required=True)
    parser.add_argument("-o", "--output", required=True)

    args = parser.parse_args()

    main(args.pi, args.recomb, args.output)
