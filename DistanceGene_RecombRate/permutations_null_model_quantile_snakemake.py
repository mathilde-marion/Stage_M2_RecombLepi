import pandas as pd
import numpy as np
import os
import argparse
from scipy import stats

# ────────────────────────────────────────────────
# Arguments
# ────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--genes", required=True)
parser.add_argument("--recomb", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--bin_size", type=int, default=1000)
parser.add_argument("--n_perm", type=int, default=100)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

np.random.seed(args.seed)
os.makedirs(args.outdir, exist_ok=True)

def tmpfile(name):
    return os.path.join(args.outdir, name)

# ────────────────────────────────────────────────
# Chargement
# ────────────────────────────────────────────────
print("Chargement des données...")
genes = pd.read_csv(args.genes, sep="\t")
recomb = pd.read_csv(args.recomb, sep="\t")
recomb["chrom"] = recomb["chrom"].astype(str).str.replace(r"b'|'", "", regex=True).str.strip()
print(f"{len(genes)} gènes | {len(recomb)} fenêtres recomb\n")

# ────────────────────────────────────────────────
# Préparer info chromosomes (via recomb)
# ────────────────────────────────────────────────
chrom_info = (
    recomb.groupby("chrom")
    .agg(chrom_min=("start", "min"), chrom_max=("end", "max"))
    .reset_index()
)
chrom_info["length"] = chrom_info["chrom_max"] - chrom_info["chrom_min"]
total_length = chrom_info["length"].sum()
chrom_info["prob"] = chrom_info["length"] / total_length

# ────────────────────────────────────────────────
# Fonctions
# ────────────────────────────────────────────────
def classify_position(pos, genes_chr):
    if len(genes_chr) == 0:
        return "intergenic_5prime", np.nan
    inside = genes_chr[(genes_chr["start"] <= pos) & (pos <= genes_chr["end"])]
    if len(inside) > 0:
        gene = inside.iloc[0]
        d1, d2 = pos - gene["start"], gene["end"] - pos
        if gene["strand"] == "+":
            return ("intragenic_5prime", d1) if d1 <= d2 else ("intragenic_3prime", d2)
        else:
            return ("intragenic_5prime", d2) if d2 <= d1 else ("intragenic_3prime", d1)
    dists = np.minimum(np.abs(pos - genes_chr["start"]), np.abs(pos - genes_chr["end"]))
    gene = genes_chr.loc[dists.idxmin()]
    d1, d2 = abs(pos - gene["start"]), abs(pos - gene["end"])
    if gene["strand"] == "+":
        return ("intergenic_5prime", d1) if d1 <= d2 else ("intergenic_3prime", d2)
    else:
        return ("intergenic_5prime", d2) if d2 <= d1 else ("intergenic_3prime", d1)

def get_rate(pos, recomb_chr):
    inside = recomb_chr[(recomb_chr["start"] <= pos) & (pos <= recomb_chr["end"])]
    if len(inside) > 0:
        return inside.iloc[0]["recombRate"], inside.iloc[0]["recomb_quantile_100"]
    before = recomb_chr[recomb_chr["end"] < pos]
    after  = recomb_chr[recomb_chr["start"] > pos]
    if len(before) == 0 or len(after) == 0:
        return np.nan, np.nan
    b, a = before.iloc[-1], after.iloc[0]
    cb, ca = (b["start"] + b["end"]) / 2, (a["start"] + a["end"]) / 2
    wb, wa = 1 / abs(pos - cb), 1 / abs(pos - ca)
    rate = (wb * b["recombRate"] + wa * a["recombRate"]) / (wb + wa)
    # Pour le quantile : on prend celui de la fenêtre la plus proche (pas d'interpolation)
    quantile = b["recomb_quantile_100"] if abs(pos - cb) <= abs(pos - ca) else a["recomb_quantile_100"]
    return rate, quantile

def run_pipeline(genes_df, recomb_df, bin_size, label=""):
    results = []
    for chrom in recomb_df["chrom"].unique():
        genes_chr  = genes_df[genes_df["chrom"] == chrom].copy().reset_index(drop=True)
        recomb_chr = recomb_df[recomb_df["chrom"] == chrom].copy().reset_index(drop=True)
        if len(recomb_chr) == 0:
            continue
        pos_min = (recomb_chr["start"].min() // bin_size) * bin_size
        pos_max = recomb_chr["end"].max()
        windows = np.arange(pos_min, pos_max + bin_size, bin_size)
        centers = windows + bin_size // 2
        classifications = [classify_position(p, genes_chr) for p in centers]
        rates_quantiles = [get_rate(p, recomb_chr) for p in centers]
        df = pd.DataFrame({
            "chrom": chrom,
            "pos": windows,
            "category": [c[0] for c in classifications],
            "distance_to_boundary": [c[1] for c in classifications],
            "recombRate": [rq[0] for rq in rates_quantiles],
            "recomb_quantile_100": [rq[1] for rq in rates_quantiles],
        })
        results.append(df)
    if not results:
        return pd.DataFrame()
    df = pd.concat(results, ignore_index=True)
    df = df[df["recombRate"].notna()]

    def refine(row):
        d = row["distance_to_boundary"]
        if pd.isna(d):
            return row["category"] + "_NA"
        s = int(d // bin_size) * bin_size
        return f"{row['category']}_{s//1000}-{(s+bin_size)//1000}kb"

    df["category_refined"] = df.apply(refine, axis=1)

    res = []
    for cat in sorted(df["category_refined"].unique()):
        sub = df[df["category_refined"] == cat]
        if len(sub) < 2:
            continue
        rates = sub["recombRate"]
        mean  = rates.mean()
        sem   = stats.sem(rates)
        ci    = stats.t.interval(0.95, len(rates)-1, loc=mean, scale=sem)
        base  = cat.rsplit("_", 1)[0]
        res.append({
            "category_refined":    cat,
            "category_base":       base,
#            "mean_rate":           mean,
            "mean_quantile_100":   sub["recomb_quantile_100"].mean(),  # <-- ajout quantile
            "ci_low":              ci[0],
            "ci_high":             ci[1],
            "mean_distance":       sub["distance_to_boundary"].mean(),
            "n_windows":           len(sub),
            "permutation":         label,
        })
    return pd.DataFrame(res)

# ────────────────────────────────────────────────
# PERMUTATIONS FULL RANDOM GENOME
# ────────────────────────────────────────────────
all_results = []
lengths = genes["end"] - genes["start"]

for perm in range(args.n_perm):
    print(f"Permutation {perm+1}/{args.n_perm}...")
    genes_perm = genes.copy()
    sampled_chroms = np.random.choice(
        chrom_info["chrom"], size=len(genes_perm), p=chrom_info["prob"]
    )
    new_starts = []
    for i, chrom in enumerate(sampled_chroms):
        row = chrom_info[chrom_info["chrom"] == chrom].iloc[0]
        chrom_min = row["chrom_min"]
        chrom_max = row["chrom_max"]
        length = lengths.iloc[i]
        start = np.random.randint(
            chrom_min, max(chrom_max - length, chrom_min + 1)
        )
        new_starts.append(start)
    genes_perm["chrom"] = sampled_chroms
    genes_perm["start"] = new_starts
    genes_perm["end"]   = genes_perm["start"] + lengths.values
    df_perm = run_pipeline(genes_perm, recomb, args.bin_size, label=f"perm_{perm+1}")
    if not df_perm.empty:
        all_results.append(df_perm)

# ────────────────────────────────────────────────
# Sauvegarde
# ────────────────────────────────────────────────
if all_results:
    out = pd.concat(all_results, ignore_index=True)
    out_path = tmpfile("permutation_null_model_quantile.tsv")
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\n✔ Sauvegardé : {out_path}")
    print(f"{len(all_results)} permutations | {len(out)} lignes")
else:
    print("⚠ Aucun résultat produit")

print("\nTerminé !")
