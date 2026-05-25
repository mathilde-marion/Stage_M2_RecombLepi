#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import argparse

# ────────────────────────────────────────────────
# Arguments Snakemake
# ────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--gff", required=True)
parser.add_argument("--recomb", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--bin_size", type=int, default=10000)
args = parser.parse_args()

path_gff = args.gff
path_recomb = args.recomb
path_out = args.outdir

os.makedirs(path_out, exist_ok=True)

def tmpfile(name):
    return os.path.join(path_out, name)

def load_or_compute(path, compute_fn):
    if os.path.exists(path):
        print(f"  → Chargement : {os.path.basename(path)}")
        return pd.read_csv(path, sep="\t")
    df = compute_fn()
    df.to_csv(path, sep="\t", index=False)
    print(f"  → Sauvegardé : {os.path.basename(path)}")
    return df

# ────────────────────────────────────────────────
# STEP 1 : GFF
# ────────────────────────────────────────────────
def compute_genes_raw():
    import gzip
    rows = []
    with gzip.open(path_gff, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            fields = line.strip().split("\t")
            if len(fields) < 9: continue
            row = fields[:8] + ["\t".join(fields[8:])]
            rows.append(row)

    gff = pd.DataFrame(rows, columns=["chrom","source","feature","start","end","score","strand","phase","attributes"])
    gff["chrom"] = gff["chrom"].astype(str)
    gff["start"] = pd.to_numeric(gff["start"])
    gff["end"]   = pd.to_numeric(gff["end"])
    return gff[gff["feature"] == "gene"].copy()

genes = load_or_compute(tmpfile("tmp_01_genes.tsv"), compute_genes_raw)

# ────────────────────────────────────────────────
# STEP 2 : recomb
# ────────────────────────────────────────────────
def compute_recomb():
    df = pd.read_csv(path_recomb, sep="\t")
    
    df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].str.replace(r"b'|'", "", regex=True)
    df["chrom"] = df["chrom"].str.strip()
    
    df["pos"] = (df["start"] + df["end"]) // 2
    return df

recomb = load_or_compute(tmpfile("tmp_02_recomb.tsv"), compute_recomb)

# ────────────────────────────────────────────────
# STEP 3 : Suppression duplicats
# ────────────────────────────────────────────────
def compute_genes_no_duplicates():
    df = genes.drop_duplicates(subset=["chrom", "start", "end", "strand"])
    return df

genes = load_or_compute(tmpfile("tmp_03_genes_no_duplicates.tsv"), compute_genes_no_duplicates)

# ────────────────────────────────────────────────
# STEP 4 : Suppression overlaps
# ────────────────────────────────────────────────
def compute_genes_no_overlaps():

    def remove_opposite(df):
        df = df.reset_index(drop=True)
        to_remove = set()
        for i in range(len(df)):
            for j in range(i+1,len(df)):
                if df.loc[i,"chrom"] != df.loc[j,"chrom"]: continue
                if df.loc[i,"strand"] == df.loc[j,"strand"]: continue
                overlap = df.loc[i,"start"] <= df.loc[j,"end"] and df.loc[j,"start"] <= df.loc[i,"end"]
                if overlap: to_remove.update([i,j])
        return df.drop(index=list(to_remove)).reset_index(drop=True)

    def remove_same(df):
        df = df.reset_index(drop=True)
        to_remove = set()
        for i in range(len(df)):
            if i in to_remove: continue
            for j in range(i+1,len(df)):
                if j in to_remove: continue
                if df.loc[i,"chrom"] != df.loc[j,"chrom"]: continue
                overlap = df.loc[i,"start"] <= df.loc[j,"end"] and df.loc[j,"start"] <= df.loc[i,"end"]
                if overlap: to_remove.add(np.random.choice([i,j]))
        return df.drop(index=list(to_remove)).reset_index(drop=True)

    df = remove_opposite(genes)
    df = remove_same(df)
    return df

genes = load_or_compute(tmpfile("tmp_04_genes_final.tsv"), compute_genes_no_overlaps)

# ────────────────────────────────────────────────
# STEP 5 : fenêtres + gène associé
# ────────────────────────────────────────────────
def compute_recomb_windows():
    window_size = int(args.bin_size)
    results = []

    def get_gene_id(attr):
        gene_id = None
        fallback_id = None

        for field in attr.split(";"):
            if field.startswith("gene_id="):
                return field.replace("gene_id=", "")
            elif field.startswith("ID="):
                fallback_id = field.replace("ID=", "")
        return gene_id if gene_id is not None else fallback_id
#        return np.nan

    def classify(pos, genes_chr):
        if len(genes_chr) == 0:
            return "intergenic_5prime", np.nan, np.nan

        inside = genes_chr[(genes_chr["start"] <= pos) & (pos <= genes_chr["end"])]

        # ── INTRAGENIC ──
        if len(inside) > 0:
            g = inside.iloc[0]
            gene_id = get_gene_id(g["attributes"])

            d1 = pos - g["start"]
            d2 = g["end"] - pos

            if g["strand"] == "+":
                if d1 <= d2:
                    return "intragenic_5prime", d1, gene_id
                else:
                    return "intragenic_3prime", d2, gene_id
            else:
                if d2 <= d1:
                    return "intragenic_5prime", d2, gene_id
                else:
                    return "intragenic_3prime", d1, gene_id

        # ── INTERGENIC ──
        dists = np.minimum(
            np.abs(pos - genes_chr["start"]),
            np.abs(pos - genes_chr["end"])
        )

        idx = dists.idxmin()
        g = genes_chr.loc[idx]
        gene_id = get_gene_id(g["attributes"])

        d1 = abs(pos - g["start"])
        d2 = abs(pos - g["end"])

        if g["strand"] == "+":
            if d1 <= d2:
                return "intergenic_5prime", d1, gene_id
            else:
                return "intergenic_3prime", d2, gene_id
        else:
            if d2 <= d1:
                return "intergenic_5prime", d2, gene_id
            else:
                return "intergenic_3prime", d1, gene_id

    for chrom in genes["chrom"].unique():
        genes_chr = genes[genes["chrom"] == chrom]

        pos_min = (genes_chr["start"].min() // window_size) * window_size
        pos_max = genes_chr["end"].max()
        windows = np.arange(pos_min, pos_max + window_size, window_size)
        centers = windows + window_size // 2

        data = [classify(p, genes_chr) for p in centers]

        results.append(pd.DataFrame({
            "chrom": chrom,
            "pos": windows,
            "category": [x[0] for x in data],
            "distance_to_boundary": [x[1] for x in data],
            "nearest_gene": [x[2] for x in data]
        }))

    return pd.concat(results, ignore_index=True)

recomb_windows = load_or_compute(
    tmpfile("tmp_05_recomb_windows_classified.tsv"),
    compute_recomb_windows
)

# ────────────────────────────────────────────────
# STEP 6 : taux recomb (interpolation pondérée)
# ────────────────────────────────────────────────
def compute_recomb_rate():
    results = []

    for chrom in recomb_windows["chrom"].unique():
        print(f"     Chromosome {chrom}")

        w = recomb_windows[recomb_windows["chrom"] == chrom].copy()
        r = recomb[recomb["chrom"] == chrom].copy().reset_index(drop=True)

        print(f"     → {len(w)} fenêtres")
        print(f"     → {len(r)} intervalles recomb")

        def get_rate(pos):
            if len(r) == 0:
                return np.nan
            # Fenêtre recomb qui contient pos
            inside = r[(r["start"] <= pos) & (pos <= r["end"])]
            if len(inside) > 0:
                return inside.iloc[0]["recombRate"]
            # Sinon interpolation pondérée entre fenêtre avant et après
            before = r[r["end"] < pos]
            after  = r[r["start"] > pos]
            if len(before) == 0 or len(after) == 0:
                return np.nan
            b  = before.iloc[-1]
            a  = after.iloc[0]
            cb = (b["start"] + b["end"]) / 2
            ca = (a["start"] + a["end"]) / 2
            wb = 1 / abs(pos - cb)
            wa = 1 / abs(pos - ca)
            return (wb * b["recombRate"] + wa * a["recombRate"]) / (wb + wa)

        w["recombRate"] = w["pos"].apply(get_rate)

        n_valid = w["recombRate"].notna().sum()
        print(f"     → {n_valid}/{len(w)} fenêtres avec taux")

        results.append(w)

    final = pd.concat(results, ignore_index=True)

    print("\n     ══ RÉSUMÉ GLOBAL ══")
    print(f"     Total fenêtres : {len(final)}")
    print(f"     Avec taux      : {final['recombRate'].notna().sum()}")
    print(f"     Sans taux      : {final['recombRate'].isna().sum()}")

    return final

recomb_final = load_or_compute(
    tmpfile("tmp_06_recomb_windows_with_rate.tsv"),
    compute_recomb_rate
)

## ────────────────────────────────────────────────
## STEP 6 : taux recomb (version robuste)
## ────────────────────────────────────────────────
#def compute_recomb_rate():
#    results = []
#
#    for chrom in recomb_windows["chrom"].unique():
#        print(f"     Chromosome {chrom}")
#
#        w = recomb_windows[recomb_windows["chrom"] == chrom].copy()
#        r = recomb[recomb["chrom"] == chrom].copy().reset_index(drop=True)
#
#        print(f"     → {len(w)} fenêtres")
#        print(f"     → {len(r)} intervalles recomb")
#
#        def get_rate(pos):
#            # Si pas de données recomb pour ce chromosome
#            if len(r) == 0:
#                return np.nan
#
#            # Distance au centre des fenêtres recomb
#            distances = np.abs(r["pos"] - pos)
#
#            # Index de la fenêtre la plus proche
#            idx = distances.idxmin()
#
#            return r.loc[idx, "recombRate"]
#
#        # Associer le taux à chaque fenêtre
#        w["recombRate"] = w["pos"].apply(get_rate)
#
#        # Petit check debug
#        n_valid = w["recombRate"].notna().sum()
#        print(f"     → {n_valid}/{len(w)} fenêtres avec taux")
#
#        results.append(w)
#
#    final = pd.concat(results, ignore_index=True)
#
#    print("\n     ══ RÉSUMÉ GLOBAL ══")
#    print(f"     Total fenêtres : {len(final)}")
#    print(f"     Avec taux      : {final['recombRate'].notna().sum()}")
#    print(f"     Sans taux      : {final['recombRate'].isna().sum()}")
#
#    return final
#
#recomb_final = load_or_compute(
#    tmpfile("tmp_06_recomb_windows_with_rate.tsv"),
#    compute_recomb_rate
#)
#
# ────────────────────────────────────────────────
# STEP 7 : catégories affinées
# ────────────────────────────────────────────────
def compute_refined():
    data = recomb_final.copy()
    bin_size = int(args.bin_size)

    def refine(row):
        d = row["distance_to_boundary"]
        if pd.isna(d): return row["category"]+"_NA"
        start = int(d//bin_size)*bin_size
        end = start+bin_size
        return f"{row['category']}_{start//1000}-{end//1000}kb"

    data["category_refined"] = data.apply(refine,axis=1)
    return data

recomb_refined = load_or_compute(tmpfile("tmp_07_recomb_windows_refined.tsv"), compute_refined)

# ────────────────────────────────────────────────
# STEP 8 : stats simples
# ────────────────────────────────────────────────
from scipy import stats

def compute_stats_basic():
    data = recomb_final[recomb_final["recombRate"].notna()]
    cats = ["intergenic_5prime","intragenic_5prime","intragenic_3prime","intergenic_3prime"]
    res = []
    for c in cats:
        sub = data[data["category"]==c]
        if len(sub)==0: continue
        rates = sub["recombRate"]
        mean = rates.mean()
        sem = stats.sem(rates)
        ci = stats.t.interval(0.95,len(rates)-1,loc=mean,scale=sem)
        res.append({
            "category": c,
            "mean_rate": mean,
            "ci_low": ci[0],
            "ci_high": ci[1],
            "sem": sem,
            "mean_distance": sub["distance_to_boundary"].mean(),
            "n_windows": len(sub)
        })
    df = pd.DataFrame(res)
    df.to_csv(os.path.join(path_out,"stats_by_category.tsv"),sep="\t",index=False)
    return df

stats_basic = compute_stats_basic()

# ────────────────────────────────────────────────
# STEP 9 : stats raffinées
# ────────────────────────────────────────────────
def compute_stats_refined():
    data = recomb_refined[recomb_refined["recombRate"].notna()]
    res = []
    for cat in sorted(data["category_refined"].unique()):
        sub = data[data["category_refined"]==cat]
        if len(sub)<2: continue
        rates = sub["recombRate"]
        mean = rates.mean()
        sem = stats.sem(rates)
        ci = stats.t.interval(0.95,len(rates)-1,loc=mean,scale=sem)
        base_category = cat.rsplit("_", 1)[0]

        res.append({
            "category_refined": cat,
            "category_base": base_category,
            "mean_rate": mean,
            "ci_low": ci[0],
            "ci_high": ci[1],
            "sem": sem,
            "mean_distance": sub["distance_to_boundary"].mean(),
            "n_windows": len(sub)
        })
    return pd.DataFrame(res)

stats_refined = compute_stats_refined()

# ────────────────────────────────────────────────
# OUTPUT FINAL (Snakemake)
# ────────────────────────────────────────────────
final_output = os.path.join(path_out,"stats_by_category_refined.tsv")
stats_refined.to_csv(final_output,sep="\t",index=False)

print(f"\nDone → {final_output}")
