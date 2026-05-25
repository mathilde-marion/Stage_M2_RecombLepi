import pandas as pd
import sys
import os

# ── Arguments ─────────────────────────────────────────
if len(sys.argv) < 8:
    print("Usage: python synteny_step4_snakemake.py pairs_list.txt recomb_list.txt output_dir p1 pop1 p2 pop2")
    sys.exit(1)

pairs_list  = sys.argv[1]
recomb_list = sys.argv[2]
out_dir     = sys.argv[3]
p1          = sys.argv[4]
pop1        = sys.argv[5]
p2          = sys.argv[6]
pop2        = sys.argv[7]

os.makedirs(out_dir, exist_ok=True)

# ── Charger listes ────────────────────────────────────
def load_pairs_list(file):
    """species -> path"""
    d = {}
    with open(file) as f:
        for line in f:
            sp, path = line.strip().split("\t")
            d[sp] = path
    return d

def load_recomb_list(file):
    """(species, pop) -> path"""
    d = {}
    with open(file) as f:
        for line in f:
            sp, pop, path = line.strip().split("\t")
            d[(sp, pop)] = path
    return d

pairs_files        = load_pairs_list(pairs_list)
recomb_files_paths = load_recomb_list(recomb_list)

# ── Charger les paires synténiques pour p1 et p2 ─────
def load_species_pairs(sp, path):
    df = pd.read_csv(path, sep="\t")
    df["pair_id"] = df["RefGene1"] + "__" + df["RefGene2"]
#    df = df.drop_duplicates(subset="pair_id").reset_index(drop=True).copy() ## AJOUT
    df = df.drop_duplicates(subset="pair_id").reset_index(drop=True).copy()
    return df

df1 = load_species_pairs(p1, pairs_files[p1])
df2 = load_species_pairs(p2, pairs_files[p2])

# ── Charger les fichiers de recomb pour (p1,pop1) et (p2,pop2) ──
def load_recomb(sp, pop, path):
    rdf = pd.read_csv(path, sep="\t")
    rdf["chrom"] = rdf["chrom"].astype(str).str.replace(r"b'|'", "", regex=True).str.strip()
    return rdf

rdf1 = load_recomb(p1, pop1, recomb_files_paths[(p1, pop1)])
rdf2 = load_recomb(p2, pop2, recomb_files_paths[(p2, pop2)])

# ── Fonction utilitaire ──────────────────────────────
def mean_recomb_in_segment(recomb_df, chrom, start, end):
    # Défense absolue contre les Series parasites
    if isinstance(start, pd.Series): start = start.iloc[0]
    if isinstance(end,   pd.Series): end   = end.iloc[0]
    if isinstance(chrom, pd.Series): chrom = chrom.iloc[0]
    start, end = float(start), float(end)
    s, e = min(start, end), max(start, end)
    mask = (
        (recomb_df["chrom"] == str(chrom)) &
        (recomb_df["start"] < e) &
        (recomb_df["end"]   > s)
    )
#    overlap = recomb_df.loc[mask, "recombRate"]
#    return overlap.mean() if len(overlap) > 0 else float("nan")
    overlap = recomb_df.loc[mask]
    if len(overlap) == 0:
        return float("nan"), float("nan")
    return overlap["recombRate"].mean(), overlap["recomb_quantile_100"].mean()

# ── Merge sur pair_id ────────────────────────────────
merged = df1.merge(df2, on="pair_id", suffixes=("_sp1", "_sp2")).reset_index(drop=True)

if merged.empty:
    print(f"⚠ Aucune paire commune entre {p1}.{pop1} et {p2}.{pop2}")
    # Snakemake attend le fichier output : on le crée vide
    pd.DataFrame().to_csv(
        os.path.join(out_dir, f"{p1}.{pop1}_VS_{p2}.{pop2}.pairs.tsv"),
        sep="\t", index=False
    )
    sys.exit(0)

# ── Sélection et renommage des colonnes ──────────────
out = merged[[
    "pair_id", "RefGene1_sp1", "RefGene2_sp1",
    "ChromosomeEspeceMap_sp1", "Gene1EspeceMap_sp1", "Gene2EspeceMap_sp1",
    "StartPosGene1EspeceMap_sp1", "EndPosGene2EspeceMap_sp1",
    "ChromosomeEspeceMap_sp2", "Gene1EspeceMap_sp2", "Gene2EspeceMap_sp2",
    "StartPosGene1EspeceMap_sp2", "EndPosGene2EspeceMap_sp2",
]].copy().rename(columns={
    "RefGene1_sp1"               : "RefGene1",
    "RefGene2_sp1"               : "RefGene2",
    "ChromosomeEspeceMap_sp1"    : f"ChromosomeEspeceMap_{p1}",
    "Gene1EspeceMap_sp1"         : f"Gene1_{p1}",
    "Gene2EspeceMap_sp1"         : f"Gene2_{p1}",
    "StartPosGene1EspeceMap_sp1" : f"StartPosGene1EspeceMap_{p1}",
    "EndPosGene2EspeceMap_sp1"   : f"EndPosGene2EspeceMap_{p1}",
    "ChromosomeEspeceMap_sp2"    : f"ChromosomeEspeceMap_{p2}",
    "Gene1EspeceMap_sp2"         : f"Gene1_{p2}",
    "Gene2EspeceMap_sp2"         : f"Gene2_{p2}",
    "StartPosGene1EspeceMap_sp2" : f"StartPosGene1EspeceMap_{p2}",
    "EndPosGene2EspeceMap_sp2"   : f"EndPosGene2EspeceMap_{p2}",
})

# ── Calcul du recombRate moyen par paire ─────────────
#for sp, pop, rdf in [(p1, pop1, rdf1), (p2, pop2, rdf2)]:
#    col_chrom = f"ChromosomeEspeceMap_{sp}"
#    col_start = f"StartPosGene1EspeceMap_{sp}"
#    col_end   = f"EndPosGene2EspeceMap_{sp}"
#    out[f"recombRate_{sp}_{pop}"] = out.apply(
#        lambda row, rdf=rdf, cc=col_chrom, cs=col_start, ce=col_end: mean_recomb_in_segment(
#            rdf, row[cc], row[cs], row[ce]
#        ), axis=1
#    )

for sp, pop, rdf in [(p1, pop1, rdf1), (p2, pop2, rdf2)]:
    col_chrom = f"ChromosomeEspeceMap_{sp}"
    col_start = f"StartPosGene1EspeceMap_{sp}"
    col_end   = f"EndPosGene2EspeceMap_{sp}"

    results = out.apply(
        lambda row, rdf=rdf, cc=col_chrom, cs=col_start, ce=col_end:
            mean_recomb_in_segment(rdf, row[cc], row[cs], row[ce]),
        axis=1
    )
    out[f"recombRate_{sp}_{pop}"]         = results.apply(lambda x: x[0])
    out[f"recomb_quantile_{sp}_{pop}"]    = results.apply(lambda x: x[1])

# ── Sauvegarde ───────────────────────────────────────
out_file = os.path.join(out_dir, f"{p1}.{pop1}_VS_{p2}.{pop2}.pairs.tsv")
out.to_csv(out_file, sep="\t", index=False)
print(f"✔ {p1}.{pop1} VS {p2}.{pop2} sauvegardé ({len(out)} paires)")
