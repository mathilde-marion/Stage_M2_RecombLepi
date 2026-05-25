## COMMENT UTILISER CE SCRIPT :  python SlidingWindowRecomb.py -i fichier_output_sortie_relernn_carte_recomb.txt -o test -s 100000 -w 1000000
## -w taille de la fenêtre glissante
## -s taille du slide

import pandas as pd
import argparse

def compute_fractional_mean_recomb_rate(input_file, output_file, window_size, slide):
    # Read the input data
    data = pd.read_csv(input_file, sep=r'\s+')
    
    results = []

    # Process each chromosome separately
    for chrom in data['chrom'].unique():
        chrom_data = data[data['chrom'] == chrom].copy()
        max_position = chrom_data['end'].max()
        #start = chrom_data['start'].min()
        start = 0 

        while start < max_position:
            end = start + window_size
            window_data = chrom_data[(chrom_data['start'] < end) & (chrom_data['end'] > start)]
            
            if not window_data.empty:
                total_recomb_rate = 0
                total_length = 0
                
                for index, row in window_data.iterrows(): #Used to deal with recombination block overlapping with the window (at the end and start)
                    block_start = max(start, row['start'])
                    block_end = min(end, row['end'])
                    block_length = block_end - block_start
                    prop_block_length=block_length/window_size
                    prop_recomb_rate = prop_block_length * row['recombRate']
                    total_length += block_length
                    total_recomb_rate += prop_recomb_rate
                    

                recomb_rate = total_recomb_rate * (window_size/total_length) if total_length > 0 else 0
                results.append({
                    'Chr': chrom,
                    'Window_Start': start,
                    'Window_End': end,
                    'Mean_recomb_rate': recomb_rate 
                })

            start += slide 
    
    # Convert results to DataFrame and save to output file
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False, sep='\t')

def main():
    parser = argparse.ArgumentParser(description="Compute mean recombination rate in sliding windows.")
    parser.add_argument('-i', required=True, help="Input file path")
    parser.add_argument('-o', required=True, help="Output file path")
    parser.add_argument('-w', type=int, required=True, help="Window size in base pairs")
    parser.add_argument('-s', type=int, required=True, help="Slide size in base pairs")

    args = parser.parse_args()

    compute_fractional_mean_recomb_rate(args.i, args.o, args.w, args.s)

if __name__ == "__main__":
    main()

print('\n\nCongrats, you have done your first sliding window analysis \\o/ !!!!!!!!!!!\n\n')
