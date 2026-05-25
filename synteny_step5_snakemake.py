import pandas as pd
import sys
import re
import numpy as np

# ── Inputs ───────────────────────────────────────────
files = sys.argv[1:-1]   # tous les fichiers en entrée
out_file = sys.argv[-1]  # dernier argument = output

# ── Fusion ───────────────────────────────────────────
merged = None

for f in files:
    df = pd.read_csv(f, sep="\t")

    if merged is None:
        merged = df
    else:
        common = ["pair_id", "RefGene1", "RefGene2"]
        new_cols = [c for c in df.columns if c not in merged.columns]

        merged = merged.merge(
            df[common + new_cols],
            on=common,
            how="outer"
        )

# ── Calcul variance par ligne ─────────────────────────
recomb_cols = [c for c in merged.columns if c.startswith("recombRate_")]
quantile_cols = [c for c in merged.columns if c.startswith("recomb_quantile_")]

# taux bruts
rates = merged[recomb_cols].replace("NA", np.nan).astype(float).replace(0, np.nan)
merged["mean_recombRate"] = rates.mean(axis=1)
merged["std_recombRate"]  = rates.std(axis=1)
merged["var_recombRate"]  = rates.var(axis=1)
merged["cv_recombRate"]   = merged["std_recombRate"] / merged["mean_recombRate"]

# taux relatifs (quantiles)
quantiles = merged[quantile_cols].replace("NA", np.nan).astype(float)
merged["mean_recomb_quantile"] = quantiles.mean(axis=1)
merged["std_recomb_quantile"]  = quantiles.std(axis=1)
merged["var_recomb_quantile"]  = quantiles.var(axis=1)
# log10
#log_rates = np.log10(rates)
#merged["mean_log_recombRate"] = log_rates.mean(axis=1)
#merged["std_log_recombRate"]  = log_rates.std(axis=1)
#merged["cv_log_recombRate"]   = merged["std_log_recombRate"] / merged["mean_log_recombRate"].abs()
#merged["n_valid_recombRate"]  = log_rates.notna().sum(axis=1)

# ── Tri naturel ──────────────────────────────────────
def natural_sort_key(s):
    return [int(x) if x.isdigit() else x for x in re.split(r'(\d+)', str(s))]

merged = merged.iloc[
    sorted(range(len(merged)),
           key=lambda i: natural_sort_key(merged["pair_id"].iloc[i]))
].reset_index(drop=True)

# ── Nettoyage ────────────────────────────────────────
merged = merged.fillna("NA")

#for col in merged.columns:
#    if not col.startswith("recombRate_") and col != "pair_id":
#        merged[col] = merged[col].apply(
#            lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() else x
#        )

for col in merged.columns:
    if not col.startswith("recombRate_") and not col.startswith("recomb_quantile_") and col != "pair_id":
        merged[col] = merged[col].apply(
            lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() else x
        )

# ── Export ───────────────────────────────────────────
merged.to_csv(out_file, sep="\t", index=False)

print(f"✔ Fichier final : {out_file}")
