"""
Integrate multiple background kmer files and generate background kmer file
"""

import gzip
import logging

def mode_makebg(arg):
    logging.info(">>>>>>>>> make_bg: Start")
    kmer_file_list = arg.kmer_files
    output_file_path = arg.output_file_path
    uniq_kmer_set = set()
    for kmer_file in kmer_file_list:
        logging.info("consider background file:{}".format(kmer_file))
        if kmer_file.endswith('gz'):
            f = gzip.open(kmer_file,'r')
            kmer_list = [line.decode().split('\t')[0] for line in f if not line.startswith(b'kmer')]
        else:
            f = open(kmer_file,'r')
            kmer_list = [line.split('\t')[0] for line in f if not line.startswith('kmer')]
        uniq_kmer_set = uniq_kmer_set.union(kmer_list)
        f.close()
    uniq_kmer_list = sorted(uniq_kmer_set)
    if arg.compressed:
        with gzip.open(output_file_path, 'wt') as f_out:
            f_out.writelines("%s\n" % l for l in uniq_kmer_list)
    else:
        with open(output_file_path,'w') as f_out:
            f_out.writelines("%s\n" % l for l in uniq_kmer_list)
    logging.info("generate unique background kmer file in {}".format(output_file_path))
    logging.info(">>>>>>>>> make_bg: Finish\n")

