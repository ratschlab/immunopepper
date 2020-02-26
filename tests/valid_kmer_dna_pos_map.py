from tests.test_immunopepper import check_kmer_pos_valid
import argparse
import sys

options = sys.argv[1:]
parser = argparse.ArgumentParser(prog='valid')
parser.add_argument("--new-junction-file",help="the output junction file outputted by filter infer dna pos."
                                               " Accept compressed and uncompressed file ", required=True)
parser.add_argument("--sample",  help="sample name to valid", required=True, default='')
parser.add_argument("--output-dir", help="specify the parent directory where mutation dict is stored if there is a pickled mutation dict",default='')
parser.add_argument("--ref-path", help="absolute path of reference genome file", required=True)
parser.add_argument("--mutation-mode", help="mutation mode (options: ref, somatic, germline, somatic_germline) [ref]",
                      required=True, default='ref')
parser.add_argument("--germline", help="absolute path of germline mutation file", required=False, default='')
parser.add_argument("--somatic", help="absolute path of somatic mutation file", required=False, default='')
args = parser.parse_args(options)

# mutation_mode='somatic_and_germline'
# new_junction_file = '{}_junction_kmer_infer.dna_pos.txt.gz'.format(mutation_mode)
# genome_file = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/sequence/genome.fa'
# sample='TCGA-A2-A0SW'
# germline_file_path='/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/germline_variants/mergedfiles_clean_stringentfilter.h5'
# somatic_file_path='/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/somatic_variants/pancan.merged.v0.2.6.PUBLIC.maf'
# '/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/immunopepper_tcga_a8f944a57a'
basic_args = ['build',
              '--samples',args.sample,
              '--splice-path','this_splicegraph',
              '--output-dir',args.output_dir,
              '--ann-path','this_ann_path',
              '--ref-path','this_ref_path',
              '--use-mut-pickle',
              ]
print(">>>> start check if {} is right".format(args.new_junction_file))
check_kmer_pos_valid(new_junction_file=args.new_junction_file,genome_file=args.ref_path,mutation_mode=args.mutation_mode,sample=args.sample,
                     germline_file_path=args.germline,somatic_file_path=args.somatic,basic_args=basic_args)
print(">>>> pass the check. All dna pos is right.")

