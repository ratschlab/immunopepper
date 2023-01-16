"""
Remove annotatation kmers
- Integrate multiple background kmer files (eg. somtic, somatic_and_germline, germline, ref annotation
- Add a column in graph (foreground) kmer file to indicate if the kmer also appear in given annotation kmer set
"""

import pandas as pd
import logging
import os
import sys

from immunopepper.io_ import read_pq_with_dict
from immunopepper.io_ import save_pd_toparquet



def mode_samplespecif(arg):
    logging.info(">>>>>>>>> Remove annotation: Start")
    bg = [os.path.exists(kmer_file) for kmer_file in arg.annot_kmer_files]
    fg = [os.path.exists(junction_kmer_file) for junction_kmer_file in arg.junction_kmer_files]
    if not sum(bg):
        logging.error("None of the sample kmer files exists, consider changing --junction_kmer_files")
        sys.exit(1)
    if not sum(fg): 
        logging.error("None of the annotation kmer files exists, consider changing --annot_kmer_file")
        sys.exit(1)
        
    bg_kmer_set = set()
    if not os.path.exists(arg.bg_file_path):
        for kmer_file in arg.annot_kmer_files:
            if os.path.exists(kmer_file):
                logging.info("...consider annotation file:{}".format(kmer_file))
                f = read_pq_with_dict(kmer_file, ['kmer'])
                bg_kmer_set.update(f['kmer'])
            else: 
                logging.info("WARNING annotation file: {} does not exist".format(kmer_file))
        save_pd_toparquet(arg.bg_file_path , pd.DataFrame(bg_kmer_set, columns = ['kmer'],  dtype='str'),
                  compression=compression, verbose=True)
        logging.info("generated unique background kmer file in {} \n ".format(arg.bg_file_path))
    else:
        bg_kmer_set = set(read_pq_with_dict(arg.bg_file_path, ['kmer'])['kmer'])
        logging.info("reading  unique background kmer file in {} \n ".format(arg.bg_file_path)
                )
    for junction_kmer_file in arg.junction_kmer_files:
        if os.path.exists(junction_kmer_file):
            logging.info("...consider foreground file:{}".format(junction_kmer_file))
            kmer_df = read_pq_with_dict(junction_kmer_file, ['kmer']).to_pandas()

            uniq_ref_kmer_set = set(kmer_df['kmer']).difference(bg_kmer_set)
            bg_flag = kmer_df['kmer'].apply(lambda x: x in uniq_ref_kmer_set)
            if  arg.remove_bg:
                kmer_df = kmer_df[bg_flag]
            else:
                kmer_df['is_neo_flag'] = bg_flag
            
            output_file_path = os.path.join(arg.output_dir, junction_kmer_file.split('/')[-1].replace('.pq', '') + '_' + arg.output_suffix + '.pq')
            save_pd_toparquet( output_file_path, kmer_df,
                  compression='SNAPPY', verbose=True)
            logging.info("output bg-removed kmer file : {} \n ".format(output_file_path))
        else:
            logging.info("WARNING foreground file: {} does not exist".format(junction_kmer_file))
    

    logging.info(">>>>>>>>> Remove annotation: Finish\n")


