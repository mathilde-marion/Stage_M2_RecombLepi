import pandas as pd
import sys
import os

if len(sys.argv) < 3:
    print("Usage: python synteny_step3.py input.tsv output.tsv")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

os.makedirs(os.path.dirname(output_file), exist_ok=True)

pairs_df = pd.read_csv(input_file, sep="\t", na_values="NA")

blocs = []
bloc_id = 0
used = [False] * len(pairs_df)

i = 0
while i < len(pairs_df):
    if used[i]:
        i += 1
        continue

    bloc = [pairs_df.iloc[i]]
    used[i] = True
    current = pairs_df.iloc[i]

    j = i + 1
    while j < len(pairs_df):
        nxt = pairs_df.iloc[j]

        if (
            not used[j] and
            nxt["RefGene1"] == current["RefGene2"] and
            nxt["ChromosomeEspeceMap"] == current["ChromosomeEspeceMap"]
        ):
            bloc.append(nxt)
            used[j] = True
            current = nxt
        j += 1

    bloc_id += 1
    first = bloc[0]
    last  = bloc[-1]

    blocs.append({
        "bloc_id": bloc_id,
        "RefGene1": first["RefGene1"],
        "RefGene2": last["RefGene2"],
        "ChromosomeEspeceMap": first["ChromosomeEspeceMap"],
        "StartPosGene1EspeceMap": first["StartPosGene1EspeceMap"],
        "EndPosGene2EspeceMap": last["EndPosGene2EspeceMap"],
        "StartPosGene1Ref": first["StartPosGene1Ref"],
        "EndPosGene2Ref": last["EndPosGene2Ref"],
        "n_genes": len(bloc) + 1,
    })

    i += 1

blocs_df = pd.DataFrame(blocs)

blocs_df.to_csv(output_file, sep="\t", index=False)
