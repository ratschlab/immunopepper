"""
Add a column in kmer tsv file to indicate if the kmer also appear in given background kmer set
"""
import pandas as pd
import logging
import gzip

def immunopepper_diff(arg):
    logging.info(">>>>>>>>> diff: Start")
    junction_kmer_file = arg.junction_kmer_file
    bg_file_path = arg.bg_file_path
    output_file = arg.output_file_path
    remove_bg = arg.remove_bg

    verbose = arg.verbose
    if bg_file_path.endswith('.gz'):
        f = gzip.open(bg_file_path,'r')
        bg_kmer_set = set(f.read().decode().split('\n'))
    else:
        f = open(bg_file_path,'r')
        bg_kmer_set = set(f.read().split('\n'))

    f.close()
    kmer_df = pd.read_csv(junction_kmer_file, sep='\t')
    uniq_ref_kmer_set = set(kmer_df['kmer']).difference(bg_kmer_set)
    bg_flag = kmer_df['kmer'].apply(lambda x: x in uniq_ref_kmer_set)
    if remove_bg:
        kmer_df = kmer_df[bg_flag]
    else:
        kmer_df['is_neo_flag'] = bg_flag
    if arg.compressed:
        kmer_df.to_csv(output_file,sep='\t',index=False,compression='gzip')
    else:
        kmer_df.to_csv(output_file, sep='\t', index=False)
    if verbose:
        logging.info("output bg-removed kmer file : {}".format(output_file))
    logging.info(">>>>>>>>> diff: Finish\n")


