"""
Remove annotatation kmers
- Integrate multiple background kmer files (eg. somtic, somatic_and_germline, germline, ref annotation
- Add a column in graph (foreground) kmer file to indicate if the kmer also appear in given annotation kmer set
"""

import pandas as pd
import logging

from .io_ import read_pq_with_dict
from .io_ import save_pd_toparquet



def mode_samplespecif(arg):
    logging.info(">>>>>>>>> Remove annotation: Start")
    uniq_kmer_set = set()
    for kmer_file in arg.kmer_files:
        logging.info("consider background file:{}".format(kmer_file))
        f = read_pq_with_dict(kmer_file, ['kmer'])
        uniq_kmer_set.update(f['kmer'])

    if arg.compressed:
        compression = 'gzip'
    else:
        compression = None

    save_pd_toparquet(arg.bg_file_path , pd.DataFrame(uniq_kmer_set, columns = ['kmer'],  dtype='str'),
                  compression=compression, verbose=True)

    logging.info("generated unique background kmer file in {}".format(arg.bg_file_path))

    logging.info("consider foreground file file:{}".format(arg.junction_kmer_file))
    kmer_df = read_pq_with_dict(arg.junction_kmer_file, ['kmer']).to_pandas()
    bg_kmer_set = set(read_pq_with_dict(arg.bg_file_path, ['kmer'])['kmer'])

    uniq_ref_kmer_set = set(kmer_df['kmer']).difference(bg_kmer_set)
    bg_flag = kmer_df['kmer'].apply(lambda x: x in uniq_ref_kmer_set)
    if  arg.remove_bg:
        kmer_df = kmer_df[bg_flag]
    else:
        kmer_df['is_neo_flag'] = bg_flag


    if arg.compressed:
        compression = 'gzip'
    else:
        compression = None

    save_pd_toparquet(arg.output_file_path, kmer_df,
                  compression=compression, verbose=True)

    logging.info("output bg-removed kmer file : {}".format(arg.output_file_path))
    logging.info(">>>>>>>>> Remove annotation: Finish\n")


