"""
Add a column in kmer tsv file to indicate if the kmer also appear in given background kmer set
"""
import pandas as pd
def immunopepper_diff(arg):
    junction_kmer_file = arg.junction_kmer_file
    bg_file_path = arg.bg_file_path
    output_file = arg.output_file_path

    verbose = arg.verbose
    with open(bg_file_path,'r') as f:
        bg_kmer_set = set(f.read().split('\n'))
    kmer_df = pd.read_csv(junction_kmer_file, sep='\t')
    uniq_ref_kmer_set = set(kmer_df['kmer']).difference(bg_kmer_set)
    kmer_df['is_neo_flag'] = kmer_df['kmer'].apply(lambda x: x in uniq_ref_kmer_set)
    kmer_df.to_csv(output_file, sep='\t', index=False)
    if verbose:
        print("add is_neo_flag in ", output_file)


