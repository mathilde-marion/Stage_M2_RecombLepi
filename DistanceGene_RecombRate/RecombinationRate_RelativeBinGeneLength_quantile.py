import pandas as pd
import numpy as np
import argparse
import re

# ── ARGUMENTS ────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--genes",  required=True)
parser.add_argument("--recomb", required=True)
parser.add_argument("--out",    required=True)
parser.add_argument("--bins",   type=int, default=10)
args = parser.parse_args()

# ── 1. LOAD ───────────────────────────────────────────────────────────────────
print("Chargement...")
genes  = pd.read_csv(args.genes,  sep="\t")
recomb = pd.read_csv(args.recomb, sep="\t")

genes["chrom"]  = genes["chrom"].astype(str)
recomb["chrom"] = recomb["chrom"].astype(str).str.replace(r"b'|'", "", regex=True).str.strip()

genes["start"]  = genes["start"].astype(int)
genes["end"]    = genes["end"].astype(int)
recomb["start"] = recomb["start"].astype(int)
recomb["end"]   = recomb["end"].astype(int)

# Colonne taux quantile
RATE_COL = "recomb_quantile_100"
if RATE_COL not in recomb.columns:
    raise ValueError(f"Colonne '{RATE_COL}' absente du fichier recomb. Colonnes disponibles : {list(recomb.columns)}")

print(f"  {len(genes)} gènes | {len(recomb)} fenêtres recomb")

# ── 2. EXTRAIRE gene_id ───────────────────────────────────────────────────────
def extract_gene_id(attr):
    attr = str(attr)
    match = re.search(r"gene_id=([^;]+)", attr)
    if match:
        return match.group(1)
    match = re.search(r"ID=([^;]+)", attr)
    if match:
        return match.group(1)
    return "NA"

genes["gene_id"] = genes["attributes"].apply(extract_gene_id)

# ── 3. CALCUL PAR GÈNE ────────────────────────────────────────────────────────
results = []

for chrom in genes["chrom"].unique():
    genes_chr  = genes[genes["chrom"] == chrom]
    recomb_chr = recomb[recomb["chrom"] == chrom].copy()

    for _, gene in genes_chr.iterrows():
        g_start  = gene["start"]
        g_end    = gene["end"]
        strand   = gene["strand"]
        gene_id  = gene["gene_id"]
        gene_len = g_end - g_start

        if gene_len <= 0:
            continue

        # Fenêtres recomb qui chevauchent ce gène
        rec_overlap = recomb_chr[
            (recomb_chr["start"] < g_end) &
            (recomb_chr["end"]   > g_start)
        ]

        for i in range(args.bins):
            bin_start = g_start + i       * gene_len / args.bins
            bin_end   = g_start + (i + 1) * gene_len / args.bins

            # Overlap pondéré avec les fenêtres recomb
            weighted_sum  = 0.0
            total_overlap = 0.0

            for _, rec in rec_overlap.iterrows():
                ov = max(0, min(bin_end, rec["end"]) - max(bin_start, rec["start"]))
                if ov > 0:
                    weighted_sum  += ov * rec[RATE_COL]
                    total_overlap += ov

            mean_rate = weighted_sum / total_overlap if total_overlap > 0 else np.nan

            # Numérotation brin-dépendante (brin - : inverser l'ordre)
            bin_index = (args.bins - i) if strand == "-" else (i + 1)

            results.append({
                "gene_id"          : gene_id,
                "chrom"            : chrom,
                "strand"           : strand,
                "gene_length"      : gene_len,
                "bin"              : bin_index,
                "recomb_quantile"  : mean_rate,
            })

# ── 4. SAVE ───────────────────────────────────────────────────────────────────
df = pd.DataFrame(results)
df.to_csv(args.out, sep="\t", index=False)
print(f"Done : {df.shape[0]} lignes → {args.out}")
print(f"  Gènes avec au moins 1 bin renseigné : {df[df['recomb_quantile'].notna()]['gene_id'].nunique()}")
print(f"  NA : {df['recomb_quantile'].isna().sum()} bins ({df['recomb_quantile'].isna().mean()*100:.1f}%)")
