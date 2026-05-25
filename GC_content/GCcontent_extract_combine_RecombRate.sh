#!/bin/bash
# ==============================================================================
# recomb_gc.sh
# Calcule le taux de recombinaison et le contenu en GC par fenêtre pour une espèce donnée
# Usage : bash recomb_gc.sh -s <species> -p <pop> [-i <input_dir>] [-g <genome>] [-o <output_dir>]
# ==============================================================================

set -euo pipefail

# --- Valeurs par défaut ---
INPUT_DIR="../../OutputRelernnCluster"
GENOME_DIR="../../genomes_ref_annotations_blast"
OUTPUT_DIR="."

# --- Aide ---
usage() {
    echo "Usage: $0 -s <species> -p <pop> [-i <input_dir>] [-g <genome_dir>] [-o <output_dir>]"
    echo ""
    echo "  -s  Nom de l'espèce (ex: lys_bel)                        [obligatoire]"
    echo "  -p  Nom de la population (ex: POP1)                       [obligatoire]"
    echo "  -i  Dossier contenant les fichiers ReLERNN  (défaut: ${INPUT_DIR})"
    echo "  -g  Dossier contenant les génomes de référence (défaut: ${GENOME_DIR})"
    echo "  -o  Dossier de sortie                        (défaut: ${OUTPUT_DIR})"
    echo ""
    echo "Exemple : bash recomb_gc.sh -s lys_bel -p POP1"
    exit 1
}

# --- Parsing des arguments ---
while getopts "s:p:i:g:o:h" opt; do
    case $opt in
        s) SPECIES="$OPTARG" ;;
        p) POP="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        g) GENOME_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# --- Vérification des arguments obligatoires ---
if [[ -z "${SPECIES:-}" || -z "${POP:-}" ]]; then
    echo "ERREUR : Les arguments -s (espèce) et -p (population) sont obligatoires."
    usage
fi

# --- Chemins dérivés ---
PREFIX="${SPECIES}.var.biallelic.fmiss80.renamed.${POP}"
INPUT_FILE="${INPUT_DIR}/${PREFIX}.PREDICT.BSCORRECTED.txt"
GENOME="${GENOME_DIR}/${SPECIES}/${SPECIES}.unmasked.fa"
OUT_PREFIX="${OUTPUT_DIR}/${SPECIES}.${POP}"

# --- Vérification des dépendances ---
for cmd in sed cut bedtools paste awk; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERREUR : '$cmd' n'est pas disponible dans le PATH."
        exit 1
    fi
done

# --- Vérification des fichiers d'entrée ---
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERREUR : Fichier introuvable : $INPUT_FILE"
    exit 1
fi
if [[ ! -f "$GENOME" ]]; then
    echo "ERREUR : Génome introuvable : $GENOME"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "  Espèce   : $SPECIES"
echo "  Pop      : $POP"
echo "  Entrée   : $INPUT_FILE"
echo "  Génome   : $GENOME"
echo "  Sortie   : $OUTPUT_DIR"
echo "========================================"

# --- Étape 1 : Nettoyage des artefacts Python b'...' ---
echo "[1/5] Nettoyage des artefacts Python (b'...')..."
CLEAN_FILE="${OUT_PREFIX}.PREDICT.BSCORRECTED_clean.txt"
sed "s/b'//g; s/'//g" "$INPUT_FILE" > "$CLEAN_FILE"

# --- Étape 2 : Extraction des colonnes BED (sans header) ---
echo "[2/5] Extraction des fenêtres BED (colonnes chrom/start/end)..."
BED_FILE="${OUT_PREFIX}.windows.bed"
tail -n +2 "$CLEAN_FILE" | cut -f1,2,3 > "$BED_FILE"

# --- Étape 3 : Calcul du contenu en GC avec bedtools nuc ---
echo "[3/5] Calcul du contenu en GC (bedtools nuc)..."
GC_RAW="${OUT_PREFIX}.gc_windows.txt"
bedtools nuc -fi "$GENOME" -bed "$BED_FILE" > "$GC_RAW"

# --- Étape 4 : Nettoyage du fichier GC ---
echo "[4/5] Extraction des colonnes GC utiles (chrom/start/end/GC)..."
GC_CLEAN="${OUT_PREFIX}.gc_clean.txt"
cut -f1,2,3,5 "$GC_RAW" > "$GC_CLEAN"

# --- Étape 5 : Fusion taux de recombinaison + GC content ---
echo "[5/5] Fusion RecombRate + GC content..."
OUT_FILE="${OUT_PREFIX}.RecombRate_GCcontent.txt"
#paste "$CLEAN_FILE" "$GC_CLEAN" \
#    | awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $5, $NF}' \
#    > "$OUT_FILE"

awk 'BEGIN{OFS="\t"}
NR==FNR { k=$1 FS $2 FS $3; gc[k]=$4; next }
NR>1 {
    k=$1 FS $2 FS $3
    if (k in gc)
        print $1,$2,$3,$5,$9,gc[k]
}' "$GC_CLEAN" "$CLEAN_FILE" > "$OUT_FILE"

# Ajout du header
sed -i '1i chrom\tstart\tend\trecombRate\tGC' "$OUT_FILE"

echo ""
echo "Fichier final : $OUT_FILE"
echo "Terminé."
