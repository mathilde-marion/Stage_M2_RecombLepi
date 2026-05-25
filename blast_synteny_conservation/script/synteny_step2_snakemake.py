import pandas as pd
import sys
import os

if len(sys.argv) < 3:
    print("Usage: python synteny_step2.py input.tsv output.tsv")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

os.makedirs(os.path.dirname(output_file), exist_ok=True)

df = pd.read_csv(input_file, sep="\t", na_values="NA")

# ── séparation ─────────────────────────────────────────
df["qchr"]   = df["qseqid"].str.split("_").str[0]
df["qgenid"] = pd.to_numeric(df["qseqid"].str.split("_").str[1], errors="coerce")

df = df.dropna(subset=["qgenid"])
df["qgenid"] = df["qgenid"].astype(int)

df_sorted = df.sort_values(["qchr", "qgenid"]).reset_index(drop=True)

# ── détection paires ───────────────────────────────────
pairs = []

for i in range(len(df_sorted) - 1):
    g1 = df_sorted.iloc[i]
    g2 = df_sorted.iloc[i + 1]

    if (
        g1["qchr"] == g2["qchr"] and
        g1["sseqid"] == g2["sseqid"] and
        abs(g2["qgenid"] - g1["qgenid"]) == 1 and
        abs(g2["sgenid"] - g1["sgenid"]) == 1
    ):
        pairs.append({
            "RefGene1": g1["qseqid"],
            "RefGene2": g2["qseqid"],
            "ChromosomeEspeceMap": g1["sseqid"],
            "Gene1EspeceMap": g1["sgenid"],
            "Gene2EspeceMap": g2["sgenid"],
            "StartPosGene1EspeceMap": int(min(g1["sstart"], g1["send"])),
            "EndPosGene2EspeceMap": int(max(g2["sstart"], g2["send"])),
            "StartPosGene1Ref": int(g1["qstartnucl"]),
            "EndPosGene2Ref": int(g2["qendnucl"]),
        })

pairs_df = pd.DataFrame(pairs)

pairs_df.to_csv(output_file, sep="\t", index=False)
