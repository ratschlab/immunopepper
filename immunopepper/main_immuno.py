

# Python libraries
import sys
import argparse
import os
import logging
from datetime import datetime

from .immunopepper_build import immunopepper_build
from .immunopepper_makebg import immunopepper_makebg
from .immunopepper_diff import immunopepper_diff
from .immunopepper_filter import immunopepper_filter

def parse_arguments(argv):

    parser = argparse.ArgumentParser(prog='immunopepper')
    subparsers = parser.add_subparsers(help='Running modes', metavar='{build, make_bg, diff, filter}')
    parser_build = subparsers.add_parser('build', help='generate kmers library from a splicegraph')
    required = parser_build.add_argument_group('MANDATORY')
    required.add_argument("--samples", nargs='+', help="list of sample names to consider", required=True, default='')
    required.add_argument("--output-dir", help="output directory [default: output]", required=True, default='output')
    required.add_argument("--ann-path", help="absolute path of reference gene annotation file", required=True)
    required.add_argument("--splice-path", help="absolute path of splicegraph file", required=True)
    required.add_argument("--ref-path", help="absolute path of reference genome file", required=True)
    required.add_argument("--mutation-mode", help="mutation mode (options: ref, somatic, germline, somatic_germline) [ref]", required=True, default='ref')

    outputs = parser_build.add_argument_group('OUTPUT OPTIONS')
    outputs.add_argument("--kmer", type=int, help="specify the k for kmer output", required=False, default=0)
    outputs.add_argument("--disable-concat",help="switch off concatenation of short exons to increase speed",action="store_true",default=False)
    outputs.add_argument("--compressed",help="compress the output files",action="store_true",default=False)

    additional_file = parser_build.add_argument_group('ADDITIONAL FILES')
    additional_file.add_argument("--germline", help="absolute path of germline mutation file", required=False, default='')
    additional_file.add_argument("--somatic", help="absolute path of somatic mutation file", required=False, default='')
    additional_file.add_argument("--count-path",help="absolute path of count hdf5 file", required=False, default=None)

    general = parser_build.add_argument_group('MISCELLANEOUS')
    general.add_argument("--process-num", metavar='N', type=int, help="Only process the first N genes in the splicegraph, default: process all", required=False, default=0)
    general.add_argument("--use-mut-pickle",help="save and use pickled mutation dict without processing the original files", action="store_true", default=False)
    general.add_argument("--verbose", type=int, help="specify the output verbosity (0 - silent, 1 - verbose, 2 - debug) [1]", required=False, default=1)

    experimental = parser_build.add_argument_group('EXPERIMENTAL')
    experimental.add_argument("--filter-redundant", help="apply redundancy filter to the exon list", action="store_true", required=False, default=False)
    #specify the absolute path to expression library sizes
    experimental.add_argument("--libsize-path", nargs='?', help=argparse.SUPPRESS,required=False, default=None)
    # specify the absolute path the the gtex_junction h5 file
    experimental.add_argument("--gtex-junction-path",help=argparse.SUPPRESS, required=False, default=None)
    # output mutated peptide even it is the same as reference peptide
    experimental.add_argument("--output-silence",help=argparse.SUPPRESS, action="store_true",default=False)
    # if count expression data is provided in h5 format, specify the code for heterzygous
    experimental.add_argument("--heter-code", type=int, help=argparse.SUPPRESS, default=0)

    parser_makebg = subparsers.add_parser('make_bg', help='integrate multiple kmer files and generate the single background kmer file')
    required = parser_makebg.add_argument_group('MANDATORY')
    required.add_argument("--kmer-files", nargs='+', help="list of kmer files output by build mode", required=True, default='')
    required.add_argument("--output-dir",help='directory to store the log file',required=True)
    required.add_argument("--output-file-path", help="output file path", required=True, default='')
    general = parser_makebg.add_argument_group('MISCELLANEOUS')
    general.add_argument("--compressed",help="compress output files",action="store_true",default=False)
    general.add_argument("--verbose", type=int, help="specify the output verbosity. Level 1 records the input and output file paths.", required=False, default=1)

    parser_diff = subparsers.add_parser('diff', help='append a new column to the junction kmer txt result file indicating if the kmer is in groundtruth')
    required = parser_diff.add_argument_group('MANDATORY')
    required.add_argument("--junction-kmer-file", help="absolute path to the foreground junction file", required=True, default='')
    required.add_argument("--bg-file-path", help="absolute path to the background file", required=True, default='')
    required.add_argument("--output-dir",help='directory to store the log file',required=True)
    required.add_argument("--output-file-path", help="output file path", required=True, default='')
    required.add_argument("--remove-bg", help="choose to simply remove background rows or add a new flag column to indicate"
                                              " if the kmer exists in the background kmers", action="store_true", required=False, default=False)

    general = parser_diff.add_argument_group('MISCELLANEOUS')
    general.add_argument("--compressed",help="compress the output files",action="store_true",default=False)
    general.add_argument("--verbose", type=int, help="specify the output verbosity. Level 1 records the output file path.", required=False, default=1)

    parser_filter = subparsers.add_parser('filter', help='apply different filter rules')
    required = parser_filter.add_argument_group('MANDATORY')
    required.add_argument("--junction-kmer-tsv-path", help="the kmer tsv file", required=True, default='')
    required.add_argument("--output-dir",help='specify the directory to store the log file',required=True)
    required.add_argument("--output-file-path", help="directory to save filtered kmer file", required=True)
    required.add_argument("--cross-junction", help="only output the cross-junction kmers", action="store_true",default=False)
    required.add_argument("--seg-expr", help="only output kmers that have segment expression greater than threshold", action="store_true",default=False)
    required.add_argument("--seg-expr-thresh", type=int, help="segment expression threshold [0]", default=0)
    required.add_argument("--junc-expr", help="only output kmers that have junction expression greater than threshold", action="store_true",default=False)
    required.add_argument("--junc-expr-thresh", type=int, help="junction expression threshold [0]", default=0)

    general = parser_filter.add_argument_group('MISCELLANEOUS')
    general.add_argument("--compressed",help="compress the output files",action="store_true",default=False)
    general.add_argument("--verbose", type=int, help="specify the output verbosity. Level 1 records the input file"
                                                     "and output file path in the log. Level 2 records the detail filter mechanism "
                                                     "and threshold.", required=False, default=1)

    if len(argv) < 1:
        parser.print_help()
        sys.exit(1)

    if len(argv) < 2:
        if argv[0] == 'build':
            parser_build.print_help()
        elif argv[0] == 'make_bg':
            parser_makebg.print_help()
        elif argv[0] == 'diff':
            parser_diff.print_help()
        elif argv[0] == 'filter':
            parser_filter.print_help()
        else:
            parser.print_help()


    pargs = parser.parse_args(argv)
    return pargs

def split_mode(options):
    arg = parse_arguments(options)
    mode = options[0]
    if not os.path.isdir(arg.output_dir):
        os.makedirs(arg.output_dir)
    now = datetime.now()
    timestamp = datetime.timestamp(now)
    runlog_name = 'run_'+mode+'_'+str(timestamp)+'.log'
    log_dir = os.path.join(arg.output_dir, runlog_name)

    file_handler = logging.FileHandler(filename=log_dir)
    if arg.verbose > 0:
        stdout_handler = logging.StreamHandler(sys.stdout)
        handlers = [file_handler, stdout_handler]
    else:
        handlers = [file_handler]
    logging.basicConfig(
                        level=logging.DEBUG,
                        handlers=handlers,
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("Command line"+str(arg))
    if mode == 'build':
        immunopepper_build(arg)
    if mode == 'make_bg':
        immunopepper_makebg(arg)
    if mode == 'diff':
        immunopepper_diff(arg)
    if mode == 'filter':
        immunopepper_filter(arg)

def cmd_entry():
    options = sys.argv[1:]
    split_mode(options)


if __name__ == "__main__":
    cmd_entry()
