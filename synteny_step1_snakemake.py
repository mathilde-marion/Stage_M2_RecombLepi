import pandas as pd
import os
import sys

if len(sys.argv) < 4:
    print("Usage: python synteny_step1.py input.tsv fasta_dir output.tsv")
    sys.exit(1)

input_file = sys.argv[1]
fasta_dir  = sys.argv[2]
output_file = sys.argv[3]

# outputs intermédiaires
out_1a = output_file.replace("step1c_nucl_coords.tsv", "step1a_sorted_by_target.tsv")
out_1b = output_file.replace("step1c_nucl_coords.tsv", "step1b_sgenid.tsv")

os.makedirs(os.path.dirname(output_file), exist_ok=True)

# ── Chargement ─────────────────────────────────────────
blast = pd.read_csv(input_file, sep="\t", na_values="NA")
blast_hits = blast.dropna(subset=["sseqid"]).copy()

# ── 1a ────────────────────────────────────────────────
blast_hits["target_pos"] = blast_hits[["sstart", "send"]].min(axis=1)
blast_hits["sseqid_sort"] = pd.to_numeric(blast_hits["sseqid"], errors="coerce")

blast_sorted = blast_hits.sort_values(
    by=["sseqid_sort", "target_pos"]
).drop(columns="sseqid_sort").reset_index(drop=True)

cols = list(blast_sorted.columns)
cols.remove("target_pos")
cols.insert(cols.index("sstart"), "target_pos")
blast_sorted = blast_sorted[cols]

blast_sorted.to_csv(out_1a, sep="\t", index=False)

# ── 1b ────────────────────────────────────────────────
blast_sorted["sgenid"] = blast_sorted.groupby("sseqid").cumcount() + 1

cols = list(blast_sorted.columns)
cols.remove("sgenid")
cols.insert(cols.index("sseqid") + 1, "sgenid")
blast_sorted = blast_sorted[cols]

blast_sorted.to_csv(out_1b, sep="\t", index=False)

# ── 1c ────────────────────────────────────────────────
fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith(".fasta")]

coords = []
for fname in fasta_files:
    qseqid = fname.split(".")[0]
    with open(os.path.join(fasta_dir, fname)) as f:
        header = f.readline().strip()

    parts = header.lstrip(">").split()
    try:
        qstartnucl = int(parts[-2])
        qendnucl   = int(parts[-1])
    except:
        continue

    coords.append({
        "qseqid": qseqid,
        "qstartnucl": qstartnucl,
        "qendnucl": qendnucl
    })

coords_df = pd.DataFrame(coords)

blast_sorted = blast_sorted.merge(coords_df, on="qseqid", how="left")

cols = list(blast_sorted.columns)
for c in ["qstartnucl", "qendnucl"]:
    cols.remove(c)

idx = cols.index("qend") + 1
cols.insert(idx, "qstartnucl")
cols.insert(idx + 1, "qendnucl")

blast_sorted = blast_sorted[cols]

blast_sorted.to_csv(output_file, sep="\t", index=False)
