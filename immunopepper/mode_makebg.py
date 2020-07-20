"""
Integrate multiple background kmer files and generate background kmer file
"""

import logging
import pandas as pd

from .io_ import read_pq_with_dict
from .io_ import save_pd_toparquet


def mode_makebg(arg):
    logging.info(">>>>>>>>> make_bg: Start")
    kmer_file_list = arg.kmer_files
    output_file_path = arg.output_file_path
    uniq_kmer_set = set()
    for kmer_file in kmer_file_list:
        logging.info("consider background file:{}".format(kmer_file))
        f = read_pq_with_dict(kmer_file, ['kmer'])
        uniq_kmer_set.update(f['kmer'])

    if arg.compressed:
        compression = 'gzip'
    else:
        compression = None

    save_pd_toparquet(output_file_path, pd.DataFrame(uniq_kmer_set, columns = ['kmer'],  dtype='str'),
                  compression=compression, verbose=True)

    logging.info("generate unique background kmer file in {}".format(output_file_path))
    logging.info(">>>>>>>>> make_bg: Finish\n")


