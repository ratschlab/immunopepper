"""
Add a column in kmer tsv file to indicate if the kmer also appear in given background kmer set
"""
import pandas as pd
import logging
import gzip

from .io_ import read_pq_with_dict
from .io_ import save_pd_toparquet


def mode_diff(arg):
    logging.info(">>>>>>>>> diff: Start")
    junction_kmer_file = arg.junction_kmer_file
    bg_file_path = arg.bg_file_path
    output_file = arg.output_file_path
    remove_bg = arg.remove_bg

    logging.info("consider foreground file file:{}".format(junction_kmer_file))
    kmer_df = read_pq_with_dict(junction_kmer_file, ['kmer']).to_pandas()
    bg_kmer_set = set(read_pq_with_dict(bg_file_path, ['kmer'])['kmer'])

    uniq_ref_kmer_set = set(kmer_df['kmer']).difference(bg_kmer_set)
    bg_flag = kmer_df['kmer'].apply(lambda x: x in uniq_ref_kmer_set)
    if remove_bg:
        kmer_df = kmer_df[bg_flag]
    else:
        kmer_df['is_neo_flag'] = bg_flag


    if arg.compressed:
        compression = 'gzip'
    else:
        compression = None

    save_pd_toparquet(output_file, kmer_df,
                  compression=compression, verbose=True)

    logging.info("output bg-removed kmer file : {}".format(output_file))
    logging.info(">>>>>>>>> diff: Finish\n")


