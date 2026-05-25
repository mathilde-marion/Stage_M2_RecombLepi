import pandas as pd
import sys

inp = sys.argv[1]
out = sys.argv[2]

# lecture correcte avec header
df = pd.read_csv(inp, sep="\t")

# garder uniquement les colonnes utiles
df = df[["chrom", "start", "end", "recombRate", "recomb_quantile_100"]]

# nettoyage chromosome b'1' -> 1
df["chrom"] = df["chrom"].astype(str).str.replace(r"^b['\"]?", "", regex=True)
df["chrom"] = df["chrom"].str.replace(r"['\"]?$", "", regex=True)

# sécurité
df = df[df["start"] < df["end"]]

# export BED
df.to_csv(out, sep="\t", index=False, header=False)
