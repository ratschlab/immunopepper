"""
Integrate multiple background kmer files and generate background kmer file
"""

import logging
import pandas as pd
import pyarrow as pa

from .io_ import save_pd_toparquet

def mode_makebg(arg):
    logging.info(">>>>>>>>> make_bg: Start")
    kmer_file_list = arg.kmer_files
    output_file_path = arg.output_file_path
    uniq_kmer_set = set()
    for kmer_file in kmer_file_list:
        logging.info("consider background file:{}".format(kmer_file))
        f = pa.parquet.read_table(kmer_file).to_pandas()
        uniq_kmer_set.update(f['kmer'].values)

    uniq_kmer_set = sorted(uniq_kmer_set)
    save_pd_toparquet(output_file_path, pd.DataFrame(uniq_kmer_set, columns = ['kmer']).head(),
                      compression=None, verbose=True)

    logging.info("generate unique background kmer file in {}".format(output_file_path))
    logging.info(">>>>>>>>> make_bg: Finish\n")


