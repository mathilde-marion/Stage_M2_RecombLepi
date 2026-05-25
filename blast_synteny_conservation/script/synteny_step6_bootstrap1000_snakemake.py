import pandas as pd
import numpy as np
from itertools import combinations
import sys

# ── Arguments Snakemake ──────────────────────────────
input_file = sys.argv[1]
out_file   = sys.argv[2]

N_PERM = 1000
SEED   = 42

# ── Chargement ──────────────────────────────────────
df = pd.read_csv(input_file, sep="\t", na_values="NA")

species = [c.replace("recombRate_", "") for c in df.columns if c.startswith("recombRate_")]

print(f"Espèces : {species}")
print(f"Paires  : {len(df)}\n")

rng = np.random.default_rng(SEED)

# ── Bootstrap ───────────────────────────────────────
results = []

for sp1, sp2 in combinations(species, 2):

    col1 = f"recombRate_{sp1}"
    col2 = f"recombRate_{sp2}"

    sub = df[["pair_id", col1, col2]].dropna()
    sub = sub[(sub[col1] > 0) & (sub[col2] > 0)].reset_index(drop=True)

    print(f"{sp1} VS {sp2} : {len(sub)} paires valides")

    if len(sub) == 0:
        continue

    pool_sp2 = sub[col2].values

    for _, row in sub.iterrows():
        x = np.log10(row[col1])
        y = np.log10(row[col2])

        obs_diff = abs(x - y)

        random_vals  = rng.choice(pool_sp2, size=N_PERM, replace=True)
        random_diffs = np.abs(x - np.log10(random_vals))

#        pval = np.mean(random_diffs < obs_diff)
        median_random = np.median(random_diffs)
        success = int(obs_diff < median_random)
        
        results.append({
            "pair_id"        : row["pair_id"],
            "sp1"            : sp1,
            "sp2"            : sp2,
            "recombRate_sp1" : row[col1],
            "recombRate_sp2" : row[col2],
            "obs_diff_log"   : round(obs_diff, 6),
            "median_rand_log": round(median_random, 6),
#            "pvalue"        : round(pval, 4),
            "success"        : success,
        })

# ── Résultats ───────────────────────────────────────
results_df = pd.DataFrame(results)
#
#def bh_correction(pvals):
#    n    = len(pvals)
#    rank = pd.Series(pvals).rank(method="first")
#    return np.minimum(pvals * n / rank.values, 1.0)
#
#results_df["pvalue_adj"] = results_df.groupby(["sp1", "sp2"])["pvalue"].transform(bh_correction)
#results_df["pvalue_adj"] = results_df["pvalue_adj"].round(4)

# ── Export ──────────────────────────────────────────
results_df.to_csv(out_file, sep="\t", index=False, na_rep="NA")

print(f"\n✔ Sauvegardé : {out_file}")
