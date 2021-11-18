
# Python libraries
import argparse
import logging
import os
import sys
import cProfile

from datetime import datetime
from mhctools.cli.args import make_mhc_arg_parser

from .mode_build import mode_build
from .mode_samplespecif import mode_samplespecif
from .mode_filter import mode_filter
from .mode_crosscohort import mode_crosscohort
from .mode_cancerspecif import mode_cancerspecif


def _add_general_args(parser):
    general = parser.add_argument_group('GENERAL')
    general.add_argument("--verbose", type=int, help="specify output verbosity (0 - warn, 1 - info, 2 - debug) [1]", required=False, default=1)
    general.add_argument("--compressed", help="compress output files", action="store_true", default=True)
    general.add_argument("--parallel", type=int, help="number of threads to be used [1]", required=False, default=1)

def parse_arguments(argv):
    parser = argparse.ArgumentParser(prog='immunopepper')
    subparsers = parser.add_subparsers(help='Running modes', metavar='{build, samplespecif, filter, cancerspecif, mhcbind}')

    ### mode_build TODO Change argument groups Mode build
    parser_build = subparsers.add_parser('build', help='generate kmers library from a splicegraph')
    required = parser_build.add_argument_group('MANDATORY')
    required.add_argument("--output-samples", nargs='+', help="list of sample names to output, names must match the graph samples", required=True, default='')
    required.add_argument("--output-dir", help="output directory [default: output]", required=True, default='output')
    required.add_argument("--ann-path", help="absolute path of reference gene annotation file", required=True)
    required.add_argument("--splice-path", help="absolute path of splicegraph file", required=True)
    required.add_argument("--ref-path", help="absolute path of reference genome file", required=True)
    required.add_argument("--mutation-mode", help="mutation mode (options: ref, somatic, germline, somatic_germline) [ref]", required=True, default='ref')
    required.add_argument("--output-fasta", help="if True outputs both the sample peptide metadata and the fasta", action="store_true", required=False, default=False)
    outputs = parser_build.add_argument_group('OUTPUT OPTIONS')
    outputs.add_argument("--kmer", nargs='+', type=int, help="list which specifys the different k for kmer output", required=False, default=[])
    outputs.add_argument("--disable-concat",help="switch off concatenation of short exons to increase speed",action="store_true",default=False)
    outputs.add_argument("--use-mut-pickle", help="save and use pickled mutation dict without processing the original files", action="store_true", default=False)
    #outputs.add_argument("--peptides_tsv", help="save the peptides outputs as tsv files instead of fasta files", action="store_true", default=False)
    outputs.add_argument("--disable_process_libsize",help="sample library size generation to increase speed",action="store_true",default=False)
    additional_file = parser_build.add_argument_group('ADDITIONAL FILES')
    additional_file.add_argument("--mutation-sample", help="sample id for the somatic and germline application, ids should match the count/graphs but equivalence with the mutation files can be specified (see --sample-name-map)", required=False, default=None)
    additional_file.add_argument("--germline", help="absolute path of germline mutation file", required=False, default='')
    additional_file.add_argument("--somatic", help="absolute path of somatic mutation file", required=False, default='')
    additional_file.add_argument("--count-path",help="absolute path of count hdf5 file", required=False, default=None)
    additional_file.add_argument("--gtex-junction-path",help="absolute path of whitelist junction file, currently only support hdf5 format. Will suport tsv"
                                                             "format in the future", required=False, default=None)
    additional_file.add_argument("--sample-name-map", help="provide an naming equivalence between the count/graphs files, the germline and the somatic mutation file Format:[ no header, 2 or 3 columns]. If 2 columns [ name in count/graphs files \t name in mutations files ] If 3 columns [name in count/graphs files \t name in germline file \t name in somatic file]", required=False, default=None) 
    _add_general_args(parser_build)
    experimental = parser_build.add_argument_group('EXPERIMENTAL')
    experimental.add_argument("--cross-graph-expr",
                          help="returns edge/segment expression matrices [kmer/peptides x samples] ",
                          action="store_true", required=False, default=False)
    experimental.add_argument("--process-chr", nargs='+',  help="Only process the list of given chromosomes from the splicegraph, default: process all", required=False, default=None)
    experimental.add_argument("--complexity-cap", type=int, help="limit the processing of the foreground to genes with complexity less than the cap", required=False, default=None)
    experimental.add_argument("--process-num", metavar='N', type=int, help="Only process the first N genes in the splicegraph, default: process all", required=False, default=0)
    experimental.add_argument("--genes-interest", help="only process the genes given as input. Path to file without header", required=False, default=None)
    experimental.add_argument("--start-id", type=int, help="development feature: start processing the graph at the given id (Tmp folder numbering in parallel mode is conserved)", required=False, default=0)
    experimental.add_argument("--filter-redundant", help="apply redundancy filter to the exon list", action="store_true", required=False, default=False)
    experimental.add_argument("--skip-annotation", help='skip the generation of the annotated peptides and kmers', action="store_true", default=False)
    experimental.add_argument("--skip-tmpfiles-rm", help='skip the removal of the intermediate directories in parallel mode', action="store_true", default=False)
    #specify the absolute path to expression library sizes
    experimental.add_argument("--libsize-path", nargs='?', help=argparse.SUPPRESS,required=False, default=None)
    # output mutated peptide even it is the same as reference peptide
    experimental.add_argument("--output-silence",help=argparse.SUPPRESS, action="store_true",default=False)
    # if count expression data is provided in h5 format, specify the code for heterzygous
    experimental.add_argument("--heter-code", type=int, help=argparse.SUPPRESS, default=0)
    experimental.add_argument("--batch-size", type=int, help="batch size for parallel processing", default=1000)
    experimental.add_argument("--all-read-frames", action="store_true", required=False, default=False)

    ### mode_samplespecif
    parser_samplespecif = subparsers.add_parser('samplespecif', help='Performs removal of the annotation to make the kmer list sample specific')
    required = parser_samplespecif.add_argument_group('MANDATORY')
    required.add_argument("--annot-kmer-files", nargs='+',  help="list of absolute paths to the annotation kmer files", required=True, default='')
    required.add_argument("--output-dir",help='directory to store the log file',required=True)
    required.add_argument("--junction-kmer-files", nargs='+',  help="absolute paths to the sample kmer files", required=True, default='')
    required.add_argument("--bg-file-path", help="absolute path to the intermediate pooled annotation file", required=True, default='')
    required.add_argument("--output-suffix", help="suffix to be append to output file path", required=True, default='no-annot')
    required.add_argument("--remove-bg", help="choose to simply remove background rows or add a new flag column to indicate"
                                              " if the kmer exists in the background kmers", action="store_true", required=False, default=False)
    _add_general_args(parser_samplespecif)

    ### mode_filter
    parser_filter = subparsers.add_parser('filter', help='apply different filter rules')
    required = parser_filter.add_argument_group('MANDATORY')
    required.add_argument("--junction-kmer-tsv-path", help="the kmer tsv file", required=True, default='')
    required.add_argument("--output-dir",help='specify the directory to store the log file',required=True)
    required.add_argument("--output-file-path", help="directory to save filtered kmer file", required=True)

    optional = parser_filter.add_argument_group('OPTIONAL')
    optional.add_argument("--cross-junction", help="only output the cross-junction kmers", action="store_true",default=False)
    optional.add_argument("--seg-expr", help="only output kmers that have segment expression greater than threshold", action="store_true",default=False)
    optional.add_argument("--seg-expr-thresh", type=float, help="segment expression threshold [0]", default=0)
    optional.add_argument("--junc-expr", help="only output kmers that have junction expression greater than threshold", action="store_true",default=False)
    optional.add_argument("--junc-expr-thresh", type=float, help="junction expression threshold [0]", default=0)
    optional.add_argument("--meta-file-path",help="specify the meta data file for more filters")
    optional.add_argument('--peptide-annotated',help="filter the kmers based on whether their original kmers appear in background peptide, 0 means keeping"
                                                     "the kmers whose original peptide does not show in background peptide. 1 means the opposite")
    optional.add_argument('--junction-annotated',help="filter the kmers based on whether their corresponding junction appear in annotation file, 0 means keeping"
                                                     "the kmers whose original junction does not show in annotation file. 1 means the opposite")
    optional.add_argument('--has-stop-codon',help="filter the kmers based on whether their corresponding sequence contains stop codon, 0 means keeping"
                                                     "the kmers whose corresponding dna does not contain stop codon. 1 means the opposite")
    optional.add_argument('--is-in-junction-list',help="filter the kmers based on whether their corresponding intron is in the junction whitelist, 0 means keeping"
                                                     "the kmers whose corresponding intron id not in the junction whitelist. 1 means the opposite")
    optional.add_argument('--is-isolated',help="filter the kmers based on whether their corresponding peptide comes from single exon or not, 0 means keeping"
                                                     "the kmers whose corresponding peptide comes from exon pairs. 1 means the opposite")
    optional.add_argument("--infer-dna-pos",help="infer the exact dna positions that output the given kmer for rna-seq filter. Need meta file provided"
                                                 "otherwise no effect",action="store_true",default=False)

    _add_general_args(parser_filter)

    ### mode_crosscohort DEPRECATED
    parser_crosscohort = subparsers.add_parser('crosscohort',
                                               help='WARNING mode deprecated: integretes kmers across cancer or normal samples. The matrices can then be used for removal of normal samples')
    required = parser_crosscohort.add_argument_group('MANDATORY')
    required.add_argument("--cores",type=int, help="number of cores for spark", required=True, default='')
    required.add_argument("--mem-per-core",type=int, help="memory per core spark", required=True)
    required.add_argument("--mutation-modes", nargs='+', help="list of all mutation modes which we would like to combine", required=True, default='')
    required.add_argument("--kmer", help='kmer', required=True)
    #required.add_argument("--output-file-path", help="directory to save filtered kmer file", required=True)
    required.add_argument("--remove-bg", help="indicate whether the background has been removed from the kmer files",
                          action="store_true", required=False, default=False)
    required.add_argument("--samples", nargs='+', help="list of all samples which we would like to combine", required=True, default='')
    required.add_argument("--input-dir", help="contains all the sample subdirectories",required=True, default='')
    required.add_argument("--output-dir", help="output directory for the integrated matrix" , required=True, default='')
    required.add_argument("--compressed_inputs", help="need to be used if .gz suffix is present on files",
                          action="store_true", required=False, default=False)
    required.add_argument("--skip-filegrouping", help="if crosscohort has ben already run once, activate to skip folder reorganisation",
                          action="store_true", required=False, default=False)
    optional = parser_crosscohort.add_argument_group('OPTIONAL')
    optional.add_argument("--output-suffix", help="suffix for the integrated matrix. e.g cancer or normals" , required=False, default='')
    _add_general_args(parser_crosscohort)

    ### mode_cancerspecif
    parser_cancerspecif = subparsers.add_parser('cancerspecif',help='Performs differential filtering against a panel of normal samples')

    rp = parser_cancerspecif.add_argument_group('Run_Parameters')
    rp.add_argument("--cores",type=int, help="number of cores", required=True, default='')
    rp.add_argument("--mem-per-core",type=int, help="memory per core", required=True)
    rp.add_argument("--parallelism", type=int, help="parallelism parameter for spark JVM", required=True, default='3')
    rp.add_argument("--out-partitions", type=int, help="number of partitions to save the final tsv file, correspond to a coalesce operation", required=False, default=None)
    rp.add_argument("--scratch-dir", help="os environement variable name containing the scratch directory path ", required=False, default='')

    so = parser_cancerspecif.add_argument_group('Samples_and_Outputs')
    so.add_argument("--kmer", help='kmer', required=True)
    so.add_argument("--output-dir", help="output directory for the filtered matrix" , required=True, default='')
    so.add_argument("--output-count", help="request to write the intermediate number of kmer at each each step to the given path (risk of slow down)" , required=False, default='')
    so.add_argument("--expression-fields-c", nargs='+', help="name of segment and junction expression field in cancer file, default ['segment_expr', 'junction_expr']",required=False, default=None)
    so.add_argument("--tag-normals", help="name for the normal cohort", required=False, default='')
    so.add_argument("--whitelist-normal", help="file containg whitelist for normal samples", required=False, default=None)
    so.add_argument("--whitelist-cancer", help="file containg whitelist for cancer samples", required=False, default=None)

    libsizes = parser_cancerspecif.add_argument_group('Libsizes')
    libsizes.add_argument("--path-cancer-libsize", help="libsize file path for cancer samples", required=False, default=None)
    libsizes.add_argument("--path-normal-libsize", help="libsize file path for normal samples", required=False, default=None)
    libsizes.add_argument("--normalizer-cancer-libsize", type=float, help="Input custom rescaling factor for cancer libsize. Default: Median of libsize", required=False,
                          default=None)
    libsizes.add_argument("--normalizer-normal-libsize", type=float, help="Input custom rescaling factor for normal libsize. Default: Median of libsize", required=False,
                          default=None)

    nsf = parser_cancerspecif.add_argument_group('Normal_Statistical_Filter')
    nsf.add_argument("--statistical", help="choose between statistical filtering or hard filtering. Default hard", action="store_true", required=False, default=False)
    nsf.add_argument("--expr-high-limit-normal", type=float, help="Normal kmers with expression >= value in >= 1 sample are truly expressed. Will not be included in statistical modelling and will be substracted from cancer set",required=False, default=None)
    nsf.add_argument("--threshold-noise-normal", type=float, help="Probability threshold on accepted noise in normals (High thresholds lead to leaner cancer kmer filtering)",required=False, default=None)
    nsf.add_argument("--tissue-grp-files", nargs='*', help="Allows the statistical modelling on normal samples to be performed on different tissue groups. Specify n paths of files, each containing the list of samples in the group. No header", required=False, default=None)

    nrf = parser_cancerspecif.add_argument_group('Normal_Recurrence_Filter')
    nrf.add_argument("--cancer-support-union", help="Choose how to combine sample-support and cohort-support in cancer. Request union, default intersection", action="store_true", required=False, default=False)
    nrf.add_argument("--id-normals",  help="Specific name for the normal cohort experiment to be saved in counting file (see output-count)", required=False, default=None)
    nrf.add_argument("--path-normal-matrix-segm", nargs='+', help="segment expression integrated matrix of kmers * samples for background", required=False, default=None)
    nrf.add_argument("--path-normal-matrix-edge", nargs='+', help="edge expression integrated matrix of kmers * samples for background", required=False, default=None)
    nrf.add_argument("--n-samples-lim-normal", type=int, help="Number of normal samples in which the expression threshold should be met", required=False, default=None)
    nrf.add_argument("--cohort-expr-support-normal", type=float, help="Expression threshold on the normal cohort applied on normalized counts of (see --expr-n-limit-normal )", required=False, default=None)
    nrf.add_argument("--path-normal-kmer-list",  nargs='+', help="Provide a list of kmer as background. Format tsv: 1 file/folder allowed,'.tsv' in filenames, format parquet: one or multiple parquets with kmer in first column", required=False, default=None)

    cbp = parser_cancerspecif.add_argument_group('Cancer_Base_Param')
    cbp.add_argument("--sample-expr-support-cancer", type=float, help="Expression threshold for the cancer target applied on normalized counts")
    cbp.add_argument("--cohort-expr-support-cancer", type=float, help="Expression threshold for the cancer cohort excluding the target sample applied on normalized counts (if --expr-n-limit-cancer and path-cancer-matrix-segm resp. edge are provided)", required=False, default=None)
    cbp.add_argument("--ids-cancer-samples", nargs='+', help=" list of all cancer samples on which to apply the filtering. If --paths-cancer-samples provided they should be given in same order", required=True, default='')
    cbp.add_argument("--mut-cancer-samples", nargs='+', help=" list of mutation modes corresponding to cancer samples. If --paths-cancer-samples provided they should be given in same order", required=True, default='')

    csf = parser_cancerspecif.add_argument_group('Cancer_Sample_Filter')
    csf.add_argument("--paths-cancer-samples", nargs='+', help="file paths of all cancer samples", required=False, default='')

    crf = parser_cancerspecif.add_argument_group('Cancer_Recurrence_Filter')
    crf.add_argument("--path-cancer-matrix-segm", nargs='+', help="segment expression integrated matrix of kmers * samples for foreground", required=False, default=None)
    crf.add_argument("--path-cancer-matrix-edge", nargs='+', help="edge expression integrated matrix of kmers * samples for foreground", required=False, default=None)
    crf.add_argument("--n-samples-lim-cancer", type=int, help="Number of normal samples in which the cohort expression threshold should be met (if --cohort-expr-support-cancer and path-cancer-matrix-segm resp. edge are provided)", required=False, default=None)

    uf = parser_cancerspecif.add_argument_group('Uniprot_Filter')
    uf.add_argument("--uniprot", help="file containg uniprot k-mers. k-mers length should match the one of the cancer and normal files", required=False, default=None)
    _add_general_args(parser_cancerspecif)
 
    ### mode_mhcbind
    parser_mhcbind = subparsers.add_parser('mhcbind',help='Perform MHC binding prediction with a wrapper for MHCtools')
    parser_mhcbind.add_argument("--mhc-software-path", help="path for MHC prediction software. e.g. netmhc3, netmhc4, netmhcpan, mhcflurry", required=True, default=None)
    parser_mhcbind.add_argument("--argstring", help="MHCtool complete command line passed as a string e.g. '--mhc-predictor netmhc3 --mhc-predictor-path /dummypath' ", required=True, default='')
    parser_mhcbind.add_argument("--output-dir", help="Specify the output directory for the MHC binding prediction", required=True, default='')
    parser_mhcbind.add_argument("--partitioned-tsv", help="Take directly the output from cancerspecif mode, need to provide the folder containing a partitionned tsv, input-peptides-file for MHC tools will be generated on the fly and saved following the path given for input-peptides-file", required=False, default=None)
    parser_mhcbind.add_argument("--bind-score-method", help="Scoring method for post binding prediction filtering e.g. score,affinity,percentile_rank for netmhc", required=False, default=None)
    parser_mhcbind.add_argument("--bind-score-threshold", type=float, help="Scoring threshold (>= threshold) for post binding prediction filtering", required=False, default=None)        
    parser_mhcbind.add_argument("--less-than", help="Scoring threshold operation becomes <= threshold", action="store_true", required=False, default=False)
    _add_general_args(parser_mhcbind)
    
    if len(argv) < 1:
        parser.print_help()
        sys.exit(1)

    if len(argv) < 2:
        if argv[0] == 'build':
            parser_build.print_help()
        elif argv[0] == 'samplespecif':
            parser_samplespecif.print_help()
        elif argv[0] == 'filter':
            parser_filter.print_help()
        elif argv[0] == "crosscohort":
            parser_crosscohort.print_help()
        elif argv[0] == "cancerspecif":
            parser_cancerspecif.print_help()
        elif argv[0] == "mhcbind":
            sys.stdout.write("------------------------------ MHCBIND IMMUNOPEPPER USAGE ------------------------------ \n \n ")
            parser_mhcbind.print_help() 
            sys.stdout.write("\n------------------------------ MHCTOOLS AVAILABLE COMMAND LINE OPTIONS ------------------------------ \n \n ")
            parser_mhc = make_mhc_arg_parser(prog="mhctools",description=("Predict MHC ligands from protein sequences"))
            parser_mhc.print_help()
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

    if arg.verbose > 0:
        stdout_handler = logging.StreamHandler(sys.stdout)
        handlers = [stdout_handler]
    ### set log level
    if arg.verbose == 0:
        log_level = logging.WARNING
    elif arg.verbose == 1:
        log_level = logging.INFO
    else:
        log_level = logging.DEBUG
        
    logging.basicConfig(
                        level=log_level,
                        handlers=handlers,
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("Command line"+str(arg))
    if mode == 'build':
        mode_build(arg)
    if mode == 'samplespecif':
        mode_samplespecif(arg)
    if mode == 'filter':
        mode_filter(arg)
    if mode == "crosscohort":
        mode_crosscohort(arg)
    if mode == "cancerspecif":
        mode_cancerspecif(arg)
    if mode == "mhcbind":
        from .mode_mhcbind import mode_mhcbind #import here due to logging conflict
        mode_mhcbind(arg)

def cmd_entry():
    #pr = cProfile.Profile()
    #pr.enable()
    options = sys.argv[1:]
    split_mode(options)
    #pr.disable()
    #pr.print_stats()
    #pr.dump_stats('/cluster/work/grlab/projects/tmp_laurie/test_memory_time_mx/cProfile.stats')

if __name__ == "__main__":
    cmd_entry()
