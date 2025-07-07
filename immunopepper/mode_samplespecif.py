"""
Remove annotatation kmers
- Integrate multiple background kmer files (eg. somtic, somatic_and_germline, germline, ref annotation
- Add a column in graph (foreground) kmer file to indicate if the kmer also appear in given annotation kmer set
"""

import pandas as pd
import logging
import os
import sys

from immunopepper.io_ import open_gz_or_normal
import glob


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
                f = open_gz_or_normal(kmer_file, 'r')
                df = pd.read_csv(f, delimiter='\t')
                bg_kmer_set.update(df['kmer'])
            else:
                logging.info("WARNING annotation file: {} does not exist".format(kmer_file))
        df_annot = pd.DataFrame(bg_kmer_set, columns = ['kmer'],  dtype='str')
        df_annot.to_csv(arg.bg_file_path, sep='\t', compression='gzip', index=False)
        logging.info("generated unique background kmer file in {} \n ".format(arg.bg_file_path))
    else:
        f = open_gz_or_normal(arg.bg_file_path, 'r')
        bg_kmer_set = set(pd.read_csv(f, delimiter='\t')['kmer'])
        logging.info("reading  unique background kmer file in {} \n ".format(arg.bg_file_path))

    for junction_kmer_file in arg.junction_kmer_files:
        if os.path.exists(junction_kmer_file):
            input_list = glob.glob('{}/*part*'.format(junction_kmer_file))
            logging.info("...consider foreground file:{}".format(input_list))
            kmer_df = pd.concat(map(lambda file: pd.read_csv(open_gz_or_normal(file, 'r'), delimiter='\t'), input_list))


            uniq_ref_kmer_set = set(kmer_df['kmer']).difference(bg_kmer_set)
            bg_flag = kmer_df['kmer'].apply(lambda x: x in uniq_ref_kmer_set)
            if  arg.remove_bg:
                kmer_df = kmer_df[bg_flag]
            else:
                kmer_df['is_neo_flag'] = bg_flag

            output_file_path = os.path.join(arg.output_dir, junction_kmer_file.split('/')[-1].replace('.gz',
                                                                                                      '') + '_' + arg.output_suffix + '.gz')

            kmer_df.to_csv(output_file_path, sep='\t', compression='gzip', index=False)
            logging.info("output bg-removed kmer file : {} \n ".format(output_file_path))
        else:
            logging.info("WARNING foreground file: {} does not exist".format(junction_kmer_file))


    logging.info(">>>>>>>>> Remove annotation: Finish\n")


