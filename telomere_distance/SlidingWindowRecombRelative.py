## COMMENT UTILISER CE SCRIPT :
## python SlidingWindowQuantiles.py \
##    -i lys_bel.var.biallelic.fmiss80.renamed.POP1.quantiles.tsv \
##    -o output.tsv \
##    -w 1000000 \
##    -s 100000
##
## -w : taille de la fenêtre glissante
## -s : taille du slide
##
## Ce script calcule la moyenne pondérée de recomb_quantile_100
## dans des fenêtres glissantes, en tenant compte des chevauchements
## partiels des blocs avec les fenêtres.

import pandas as pd
import argparse

def compute_fractional_mean_quantile(input_file, output_file, window_size, slide):

    # Lecture du fichier
    data = pd.read_csv(input_file, sep=r'\s+')

    results = []

    # Analyse chromosome par chromosome
    for chrom in data['chrom'].unique():

        chrom_data = data[data['chrom'] == chrom].copy()

        max_position = chrom_data['end'].max()

        start = 0

        while start + window_size < max_position:

            end = start + window_size

            # Blocs qui chevauchent la fenêtre
            window_data = chrom_data[
                (chrom_data['start'] < end) &
                (chrom_data['end'] > start)
            ]

            if not window_data.empty:

                total_quantile = 0
                total_length = 0

                # Gestion des chevauchements partiels
                for index, row in window_data.iterrows():

                    block_start = max(start, row['start'])
                    block_end = min(end, row['end'])

                    block_length = block_end - block_start

                    # Contribution pondérée du bloc
                    prop_block_length = block_length / window_size
                    prop_quantile = (
                        prop_block_length * row['recomb_quantile_100']
                    )

                    total_length += block_length
                    total_quantile += prop_quantile

                # Normalisation si la fenêtre n'est pas entièrement couverte
                mean_quantile = (
                    total_quantile * (window_size / total_length)
                    if total_length > 0 else 0
                )

                results.append({
                    'Chr': chrom,
                    'Window_Start': start,
                    'Window_End': end,
                    'Mean_recomb_quantile_100': mean_quantile
                })

            start += slide

    # Sauvegarde
    results_df = pd.DataFrame(results)

    results_df.to_csv(
        output_file,
        index=False,
        sep='\t'
    )

def main():

    parser = argparse.ArgumentParser(
        description="Compute mean recombination quantile in sliding windows."
    )

    parser.add_argument(
        '-i',
        required=True,
        help="Input file path"
    )

    parser.add_argument(
        '-o',
        required=True,
        help="Output file path"
    )

    parser.add_argument(
        '-w',
        type=int,
        required=True,
        help="Window size in base pairs"
    )

    parser.add_argument(
        '-s',
        type=int,
        required=True,
        help="Slide size in base pairs"
    )

    args = parser.parse_args()

    compute_fractional_mean_quantile(
        args.i,
        args.o,
        args.w,
        args.s
    )

if __name__ == "__main__":
    main()

print(
    '\n\nCongrats, you have done your sliding window quantile analysis \\o/ !!!!!!!!!!!\n\n'
)
