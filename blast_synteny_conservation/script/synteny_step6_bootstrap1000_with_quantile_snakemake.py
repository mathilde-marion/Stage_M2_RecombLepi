import pandas as pd
import numpy as np
from itertools import combinations
import sys

input_file = sys.argv[1]
out_file   = sys.argv[2]
N_PERM = 1000
SEED   = 42

df = pd.read_csv(input_file, sep="\t", na_values="NA")

recomb_cols   = [c.replace("recombRate_", "")      for c in df.columns if c.startswith("recombRate_")]
quantile_cols = [c.replace("recomb_quantile_", "") for c in df.columns if c.startswith("recomb_quantile_")]

# vérification que les deux listes sont cohérentes
assert set(recomb_cols) == set(quantile_cols), "Colonnes recombRate et recomb_quantile non concordantes"
species = recomb_cols

print([c for c in df.columns if "recomb" in c.lower()])

print(f"Espèces : {species}")
print(f"Paires  : {len(df)}\n")

rng = np.random.default_rng(SEED)

# ── Fonction bootstrap générique ─────────────────────
def run_bootstrap(df, species, col_prefix, log_transform):
    results = []
    for sp1, sp2 in combinations(species, 2):
        col1 = f"{col_prefix}{sp1}"
        col2 = f"{col_prefix}{sp2}"
        sub = df[["pair_id", col1, col2]].dropna()
        sub = sub[(sub[col1] > 0) & (sub[col2] > 0)].reset_index(drop=True)
        print(f"[{col_prefix}] {sp1} VS {sp2} : {len(sub)} paires valides")
        if len(sub) == 0:
            continue
        pool_sp2 = sub[col2].values
        for _, row in sub.iterrows():
            x = np.log10(row[col1]) if log_transform else row[col1]
            y = np.log10(row[col2]) if log_transform else row[col2]
            obs_diff = abs(x - y)
            random_vals  = rng.choice(pool_sp2, size=N_PERM, replace=True)
            random_diffs = np.abs(x - (np.log10(random_vals) if log_transform else random_vals))
            median_random = np.median(random_diffs)
            success = int(obs_diff < median_random)
            results.append({
                "pair_id"        : row["pair_id"],
                "sp1"            : sp1,
                "sp2"            : sp2,
                f"{col_prefix}sp1" : row[col1],
                f"{col_prefix}sp2" : row[col2],
                "obs_diff"       : round(obs_diff, 6),
                "median_rand"    : round(median_random, 6),
                "success"        : success,
            })
    return pd.DataFrame(results)

# ── Bootstrap taux absolu (log-transformé) ───────────
results_abs = run_bootstrap(df, species, col_prefix="recombRate_", log_transform=True)
results_abs.columns = [c if c in ("pair_id","sp1","sp2") 
                        else c.replace("obs_diff", "obs_diff_log")
                             .replace("median_rand", "median_rand_log")
                        for c in results_abs.columns]

# ── Bootstrap taux relatif (quantile, pas de log) ────
results_quant = run_bootstrap(df, species, col_prefix="recomb_quantile_", log_transform=False)
results_quant.columns = [c if c in ("pair_id","sp1","sp2")
                          else c.replace("obs_diff", "obs_diff_quantile")
                               .replace("median_rand", "median_rand_quantile")
                               .replace("success", "success_quantile")
                          for c in results_quant.columns]

# ── Fusion des deux résultats ─────────────────────────
merge_cols = ["pair_id", "sp1", "sp2"]
results_df = results_abs.merge(
    results_quant.drop(columns=["recomb_quantile_sp1", "recomb_quantile_sp2"]),
    on=merge_cols,
    how="outer"
)

results_df.to_csv(out_file, sep="\t", index=False, na_rep="NA")
print(f"\n✔ Sauvegardé : {out_file}")
