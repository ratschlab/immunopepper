
# Python libraries
import argparse
import logging
import os
import sys

from datetime import datetime
from mhctools.cli.args import make_mhc_arg_parser

from immunopepper.mode_build import mode_build
from immunopepper.mode_samplespecif import mode_samplespecif
from immunopepper.mode_cancerspecif import mode_cancerspecif

def _add_general_args(parser):
    general = parser.add_argument_group('GENERAL')
    general.add_argument("--verbose", type=int, help="specify output verbosity (0 - warn, 1 - info, 2 - debug) [1]", required=False, default=1)

def parse_arguments(argv):
    parser = argparse.ArgumentParser(prog='immunopepper')
    subparsers = parser.add_subparsers(help='Running modes', metavar='{build, samplespecif, cancerspecif, mhcbind}')

    ### mode_build
    parser_build = subparsers.add_parser('build', help='generate kmers library from a splicegraph')

    required = parser_build.add_argument_group('MANDATORY')
    required.add_argument("--output-dir", help="output directory [default: output]", required=True, default='output')
    required.add_argument("--ann-path", help="absolute path of reference gene annotation file", required=True)
    required.add_argument("--splice-path", help="absolute path of splicegraph file", required=True)
    required.add_argument("--ref-path", help="absolute path of reference genome file", required=True)
    required.add_argument("--kmer", nargs='+', type=int, help="list which specifies the different k for kmer output", required=True, default=[])

    submodes = parser_build.add_argument_group('SUBMODES PROCESS: Conceptual choices about the processing required')
    submodes.add_argument("--libsize-extract", help="goes through the graph only to output gene quantifications and library sizes. Skips neoepitope generation", action="store_true", required=False, default=False)
    submodes.add_argument("--all-read-frames", help="switch from annotated reading frames to exhaustive translation", action="store_true", required=False, default=False)
    submodes.add_argument("--count-path", help="absolute path of count hdf5 file. If not provided candidates are output without expression quantification", required=False, default=None)
    submodes.add_argument("--output-samples", nargs='+', help="list of sample names to output, names must match the graph samples. If not provided, all count samples will be output (runs faster)", required=False, default=[])
    submodes.add_argument("--heter-code", type=int, help=argparse.SUPPRESS, default=0) # if count expression data is provided in h5 format, specify the code for heterzygous

    parameters = parser_build.add_argument_group('TECHNICAL PARAMETERS: Parameters to dump intermediate data or to optimize processing')
    parameters.add_argument("--compressed", help="compress output files", action="store_true", default=True)
    parameters.add_argument("--parallel", type=int, help="number of threads to be used [1]", required=False, default=1)
    parameters.add_argument("--batch-size", type=int, help="batch size for parallel processing", default=1000)
    parameters.add_argument("--pickle-samples", nargs='+', help="list of sample names to pickle, ids should match the count/graphs but equivalence with the mutation files can be specified (see --sample-name-map)", required=False, default=[])

    subset = parser_build.add_argument_group('SUBSET: Process subsets of the graph')
    subset.add_argument("--process-chr", nargs='+',  help="Only process the list of given chromosomes from the splicegraph, default: process all", required=False, default=None)
    subset.add_argument("--complexity-cap", type=int, help="limit the processing of the foreground to genes with complexity less than the cap", required=False, default=None)
    subset.add_argument("--genes-interest", help="only process the genes given as input. One gene-ID per line, no header.", required=False, default=None)
    subset.add_argument("--start-id", type=int, help="development feature: start processing the graph at the given id (Tmp folder numbering in parallel mode is conserved)", required=False, default=0)
    subset.add_argument("--process-num", metavar='N', type=int, help="Only process the first N genes in the splicegraph, default: process all", required=False, default=0)

    outputs = parser_build.add_argument_group('OUTPUT OPTIONS: Optional choices about the output formatting and filtering')
    outputs.add_argument("--skip-annotation", help='skip the generation of the annotated peptides and kmers', action="store_true", default=False)
    outputs.add_argument("--keep-tmpfiles", help='keep the intermediate directories in parallel mode', action="store_true", default=False)
    outputs.add_argument("--libsize-path", nargs='?', help=argparse.SUPPRESS, required=False, default=None)     #specify the absolute path to expression library sizes if we want to append to a file
    outputs.add_argument("--output-fasta", help="if True, outputs both the sample peptide metadata and the fasta, else outputs only metadata", action="store_true", required=False, default=False)
    outputs.add_argument("--force-ref-peptides", help="output mutated peptide even if it is the same as reference peptide", action="store_true", default=False)
    outputs.add_argument("--filter-redundant", help="apply redundancy filter to the exon list", action="store_true", required=False, default=False)
    outputs.add_argument("--kmer-database", help="absolute path of database file with kmers on one column (no header). If provided, the kmers matching the given database with I and L equivalence will not be output.", required=False, default=None)
    outputs.add_argument("--gtex-junction-path", help="absolute path of whitelist junction file, currently only support hdf5 format.", required=False, default=None)
    parameters.add_argument("--disable-concat", help="switch off concatenation of short exons to increase speed", action="store_true", default=False)
    parameters.add_argument("--disable-process-libsize", help="sub-sample library size generation to increase speed", action="store_true", default=False)

    mutation = parser_build.add_argument_group('MUTATION FILES: Arguments needed for integration of the germline variants or somatic mutations')
    mutation.add_argument("--mutation-sample", help="sample id for the somatic and germline application, ids should match the count/graphs but equivalence with the mutation files can be specified (see --sample-name-map)", required=False, default=None)
    mutation.add_argument("--germline", help="Absolute path to an optional germline VCF or MAF mutation file", required=False, default='')
    mutation.add_argument("--somatic", help="Absolute path to an optional somatic VCF or MAF mutation file", required=False, default='')
    mutation.add_argument("--sample-name-map", help="provide a name mapping between the count/graphs files, the germline and the somatic mutation file Format:[ no header, 2 or 3 columns]. If 2 columns [ name in count/graphs files \t name in mutations files ] If 3 columns [name in count/graphs files \t name in germline file \t name in somatic file]", required=False, default=None)
    mutation.add_argument("--use-mut-pickle", help="save and use pickled mutation dict without processing the original files", action="store_true", default=False)

    _add_general_args(parser_build)

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

    ### mode_cancerspecif
    parser_cancerspecif = subparsers.add_parser('cancerspecif', help='Performs differential filtering against a panel of normal samples')

    spark = parser_cancerspecif.add_argument_group('TECHNICAL PARAMETERS: parameters to control spark processing')
    spark.add_argument("--cores", type=int, help="number of cores", required=True, default='')
    spark.add_argument("--mem-per-core", type=int, help="memory per core", required=True)
    spark.add_argument("--parallelism", type=int, help="parallelism parameter for spark JVM", required=True, default='3')
    spark.add_argument("--out-partitions", type=int, help="number of partitions to save the final tsv file, correspond to a coalesce operation", required=False, default=None)
    spark.add_argument("--scratch-dir", help="os environement variable name containing the cluster scratch directory path", required=False, default='')
    spark.add_argument("--interm-dir-norm", help="custom scatch dir path to save intermediate cancer files", required=False, default='')
    spark.add_argument("--interm-dir-canc", help="custom scatch dir path to save intermediate normal files", required=False, default='')

    helpers = parser_cancerspecif.add_argument_group('INPUT HELPERS: Help the software understand the input files')
    helpers.add_argument("--kmer", help='kmer', required=True)
    helpers.add_argument("--expression-fields-c", nargs='+', help="name of segment and junction expression field in cancer file, default ['segment_expr', 'junction_expr']",required=False, default=None)
    helpers.add_argument("--ids-cancer-samples", nargs='+',help=" list of all cancer samples on which to apply the filtering. If --paths-cancer-samples provided they should be given in same order", required=True, default='')
    helpers.add_argument("--mut-cancer-samples", nargs='+', help=" list of mutation modes corresponding to cancer samples. If --paths-cancer-samples provided they should be given in same order", required=True, default='')

    inputs = parser_cancerspecif.add_argument_group('GENERAL INPUT FILES: Files and parameters to be provided to the software regardless of the filtering strategy')
    inputs.add_argument("--whitelist-normal", help="file containg whitelist for normal samples", required=False, default=None)
    inputs.add_argument("--whitelist-cancer", help="file containg whitelist for cancer samples", required=False, default=None)
    inputs.add_argument("--path-cancer-libsize", help="libsize file path for cancer samples", required=False, default=None)
    inputs.add_argument("--path-normal-libsize", help="libsize file path for normal samples", required=False, default=None)
    inputs.add_argument("--normalizer-cancer-libsize", type=float, help="Input custom rescaling factor for cancer libsize. Default: Median of libsize", required=False, default=None)
    inputs.add_argument("--normalizer-normal-libsize", type=float, help="Input custom rescaling factor for normal libsize. Default: Median of libsize", required=False, default=None)

    outputs = parser_cancerspecif.add_argument_group('GENERAL OUTPUT FILES: Files and parameters to be output by the software regardless of the filtering strategy')
    outputs.add_argument("--output-dir", help="output directory for the filtered matrix" , required=True, default='')
    outputs.add_argument("--output-count", help="request to write the intermediate number of kmer at each each step to the given path (risk of slowdown)" , required=False, default='')
    outputs.add_argument("--tag-normals", help="name for the normal cohort output files, use when various normal cohorts", required=False, default='')
    outputs.add_argument("--tag-prefix", help="prefix to use for the output files, use when several conditions", required=False, default='')

    nsf = parser_cancerspecif.add_argument_group('NORMAL SAMPLES: "Submode Statistical Filter". Fits a NB distribution on normal kmers and use a probabilistic threshold for normal background inclusion')
    nsf.add_argument("--statistical", help="choose between statistical filtering or hard filtering. Default hard", action="store_true", required=False, default=False)
    nsf.add_argument("--expr-high-limit-normal", type=float, help="Normal kmers with expression >= value in >= 1 sample are truly expressed. Will not be included in statistical modelling and will be substracted from cancer set",required=False, default=None)
    nsf.add_argument("--threshold-noise-normal", type=float, help="Probability threshold on accepted noise in normals (High thresholds lead to leaner cancer kmer filtering)",required=False, default=None)
    nsf.add_argument("--tissue-grp-files", nargs='*', help="Allows the statistical modelling on normal samples to be performed on different tissue groups. Specify n paths of files, each containing the list of samples in the group. No header", required=False, default=None)

    nrf = parser_cancerspecif.add_argument_group('NORMAL SAMPLES: "Submode Recurrence hard Filter". Normal background inclusion is based on a combination of two filters (a) Number of reads to be expressed in any sample (b) Number of samples to have any read. The filters can be requested independently.')
    nrf.add_argument("--path-normal-matrix-segm", nargs='+', help="Segment expression integrated matrix of kmers * samples for background", required=False, default=None)
    nrf.add_argument("--path-normal-matrix-edge", nargs='+', help="Edge expression integrated matrix of kmers * samples for background", required=False, default=None)
    nrf.add_argument("--n-samples-lim-normal", type=int, help="Number of normal samples in which any number of reads is required (>0)", required=False, default=None)
    nrf.add_argument("--cohort-expr-support-normal", type=float, help="Normalized expression threshold for the normal cohort required in any sample (>=1)", required=False, default=None)

    crf = parser_cancerspecif.add_argument_group('CANCER SAMPLES: "Submode Recurrence and support filter". Cancer foreground inclusion is based on a combination of two filters (a) Number of reads to be expressed in the cancer "target" sample (b) Number of reads to be expressed in more than n samples in other cancer cohort samples. The filters can be requested independently. Note that use of (b) requires matrix input files.' )
    crf.add_argument("--sample-expr-support-cancer", type=float, help="Normalized expression threshold for the cancer target sample. (Can always be specified)")
    crf.add_argument("--cohort-expr-support-cancer", type=float, help="Normalized expression threshold for the cancer cohort excluding the target sample which should be met in n samples (if --expr-n-limit-cancer and path-cancer-matrix-segm resp. edge are provided)", required=False, default=None)
    crf.add_argument("--n-samples-lim-cancer", type=int, help="Number of cancer samples in which the cancer cohort expression threshold should be met (if --cohort-expr-support-cancer and path-cancer-matrix-segm resp. edge are provided)", required=False, default=None)
    crf.add_argument("--path-cancer-matrix-segm", nargs='+', help="List of cohort cancer matrix files containing segment expression in [kmers x samples]", required=False, default=None)
    crf.add_argument("--path-cancer-matrix-edge", nargs='+', help="List of cohort cancer matrix files containing edge expression in [kmers x samples]", required=False, default=None)
    crf.add_argument("--cancer-support-union", help="Choose how to combine sample-expr and cohort-expr support in cancer. Request union, default intersection", action="store_true", required=False, default=False)

    more_backgrouds = parser_cancerspecif.add_argument_group('ADDITIONAL BACKGROUNDS: Other lists of kmers to be removed')
    more_backgrouds.add_argument("--path-normal-kmer-list",  nargs='+', help="Provide a list of kmer to be used independantly as background. Format tsv: 1 file/folder allowed,'.tsv' in filenames, format parquet: one or multiple parquets with kmer in first column", required=False, default=None)
    more_backgrouds.add_argument("--uniprot", help="file containg uniprot k-mers. k-mers length should match the one of the cancer and normal files", required=False, default=None)

    more_filters = parser_cancerspecif.add_argument_group('ADDITIONAL FILTERS: Currently filters on the annotation status of the junction coordinate and the reading frame of the kmers are supported')
    more_filters.add_argument("--filterNeojuncCoord", choices=['C', 'N', 'A'], required=False, default='', help="Retain kmers generated from neojunctions i.e. whose junction coordinates are not found in the annotation. Values: 'C', 'N', 'A' to perform filtering in cancer, normal or both sets respectively")
    more_filters.add_argument("--filterAnnotatedRF", choices=['C', 'N', 'A'], required=False, default='', help="Retain kmers generated from annotated reading frames i.e. whose reading frames are taken from annotated transcript and not propagated through the graph. Values: 'C', 'N', 'A' to perform filtering in cancer, normal or both sets respectively")


    development = parser_cancerspecif.add_argument_group('DEVELOPMENT PARAMETERS')
    development.add_argument("--tot-batches", type=int, help="Filter foreground and in background kmers based on hash function. Set number of batches",required=False, default=None)
    development.add_argument("--batch-id", type=int, help="Filter foreground and in background kmers based on hash function. Set 0<= batch_id <tot_batches",required=False, default=None)
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
