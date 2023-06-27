import glob
import logging
import os
import pandas as pd
import shlex
import sys

from mhctools.cli.script import main as mhctools_main

def mode_mhcbind(arg):
    args_list=shlex.split(arg.argstring)
    if '--output-csv' not in args_list:
        logging.error('output-csv was not provided in --arg.argstring')
        sys.exit(1)

    ###Process the outputs from cancerspecif mode
    if '--input-peptides-file' in args_list and arg.partitioned_tsv:
        logging.info("Process the outputs from cancerspecif mode")
        input_peptides_file_idx = [i+1 for i,j in enumerate(args_list) if j == '--input-peptides-file']
        input_peptides_file = args_list[input_peptides_file_idx[0]]
        input_list = glob.glob('{}/*part*'.format(arg.partitioned_tsv))
        if len(input_list) == 0:
             logging.error("No file partitions in {}. Please check --partitioned-tsv".format(arg.partitioned_tsv))
             sys.exit(1)
        input_= pd.concat(map(lambda file: pd.read_csv(file , sep ='\t'), input_list))
        input_[['kmer']].to_csv(input_peptides_file, header = None, index = None)

    ### Launch MHC Tools
    logging.info("Launch MHC Tools with command {}".format(arg.argstring))
    os.environ["PATH"] +=":{}".format(arg.mhc_software_path)
    mhctools_main(args_list)

    output_csv_idx = [i+1 for i,j in enumerate(args_list) if j == '--output-csv']
    output_csv_path = args_list[output_csv_idx[0]]
    output_csv = pd.read_csv(output_csv_path, sep = ',' )
    ### Merge the outputs from MHC bind and cancerspecif modes
    if '--input-peptides-file' in args_list and arg.partitioned_tsv:
        output_csv = input_.merge( output_csv, left_on ='kmer', right_on = 'peptide').drop('peptide', axis = 1)

    ### Perform some rank or affinity binding filtering
    if (arg.bind_score_method is not None) and (arg.bind_score_threshold is not None):
        if arg.less_than:
            logging.info("Perform filtering with {} <= {}".format(arg.bind_score_method, arg.bind_score_threshold))
            output_csv = output_csv.loc[output_csv[arg.bind_score_method] <= arg.bind_score_threshold, :]
            saving_suffix='_With{}LessLim{}'.format(arg.bind_score_method, arg.bind_score_threshold)
        else:
            logging.info("Perform filtering with {} >= {}".format(arg.bind_score_method, arg.bind_score_threshold))
            output_csv = output_csv.loc[output_csv[arg.bind_score_method] >= arg.bind_score_threshold, :]
            saving_suffix='_With{}MoreLim{}'.format(arg.bind_score_method, arg.bind_score_threshold)
    else:
        saving_suffix = ''

    if 'tsv' or 'csv' in output_csv:
        base_out = '.'.join(output_csv_path.split('.')[:-1])
    else:
        base_out = output_csv_path
    format_tag = '.tsv'

    logging.info("Saving to {}".format(base_out + saving_suffix + format_tag))
    output_csv.to_csv( base_out + saving_suffix + format_tag, sep = '\t', index = None)
