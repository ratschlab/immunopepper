from __future__ import print_function

# Python libraries
import sys
import cPickle
import os
import timeit
import argparse
import gzip

# External libraries
import Bio.SeqIO as BioIO
import h5py

# immuno module
from immuno_print import print_memory_diags
from immuno_preprocess import genes_preprocess,preprocess_ann,parse_gene_metadata_info,parse_mutation_from_maf,parse_mutation_from_vcf_h5,parse_junction_meta_info,parse_mutation_from_vcf
from immuno_mutation import get_mutation_mode_from_parser
from immuno_model import annotate_gene_opt
from immuno_filter import get_filtered_output_list
from modules.classes import gene

### Example usage
# python main_immuno.py --samples TCGA-13-1489 --output_dir t --splice_path quick_test_data/sample_gene.pkl --ann_path quick_test_data/small.gtf --ref_path quick_test_data/smallgene34.fa
# python exon_filter.py --meta_file t/TCGA-13-1489/ref_meta1.tsv.gz --peptide_file t/TCGA-13-1489/ref_peptide2.fa

def parse_arguments(argv):

    parser = argparse.ArgumentParser(argv)
    parser.add_argument("--samples", nargs='+', help="the sample names, can specify more than one sample", required=False, default='')
    parser.add_argument("--output_dir", help="specify the output directory [default: test]", required=False, default='test')
    parser.add_argument("--ann_path", help="specify the absolute path of annotation file", required=False)
    parser.add_argument("--splice_path", help="specify the absolute path of splicegraph file", required=False)
    parser.add_argument("--ref_path", help="specify the absolute path of reference gene file to the work_dir", required=False)
    parser.add_argument("--libsize_path", nargs='?', help="specify the absolute path to expression library sizes",required=False, default=None)
    parser.add_argument("--vcf_path", nargs=1, help="specify the absolute path of vcf file", required=False, default='')
    parser.add_argument("--maf_path", nargs=1, help="specify the absolute path of maf file", required=False, default='')
    parser.add_argument("--count_path",help="specify the absolute path of the count h5 file", required=False, default=None)
    parser.add_argument("--gtex_junction_path",help="specify the absolute path the the gtex_junction h5 file", required=False, default=None)
    parser.add_argument("--process_num", type=int, help="Only process the first *process_num* gene in the splicegraph,default,0, means process all", required=False, default=0)
    parser.add_argument("--is_filter", help="apply redundancy filter to the exon list", action="store_false", required=False, default=False)
    parser.add_argument("--debug", help="generate debug output", action="store_true", required=False, default=False)


    # if len(argv) < 2:
    #     parser.print_help()
    #     sys.exit(1)

    pargs = parser.parse_args()
    return pargs


def main():
    print(os.path.abspath(os.curdir))
    arg = parse_arguments(sys.argv)

    ## for debugging in pycharm
    arg.output_dir = '../tests'
    arg.ann_path = '../tests/data/test1.gtf'
    arg.ref_path = '../tests/data/test1.fa'
    arg.splice_path = '../tests/data/spladder/genes_graph_conf3.merge_graphs.pickle'
    #arg.gtex_junction_path = '../tests/data/gtex_junctions.hdf5'
    #arg.vcf_path = ['../tests/data/test1pos.vcf']
    #arg.maf_path = ['../tests/data/test1.maf']
    arg.count_path = '../tests/data/spladder/genes_graph_conf3.merge_graphs.count.hdf5'
    arg.samples = ['test1']
    mutation_mode, vcf_file_path, maf_file_path = get_mutation_mode_from_parser(arg)

    # load genome sequence data
    seq_dict = {}
    start_time = timeit.default_timer()
    interesting_chr = map(str, range(1, 23)) + ["X", "Y", "MT"]

    # read reference genome for standard chromosomes
    print('Parsing genome sequence ...')
    for record in BioIO.parse(arg.ref_path, "fasta"):
        if record.id in interesting_chr:
            seq_dict[record.id] = str(record.seq).strip()
    end_time = timeit.default_timer()
    print('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()

    # read and process the annotation file
    print('Building lookup structure ...')
    start_time = timeit.default_timer()
    gene_cds_begin_dict, gene_to_transcript_table, transcript_to_cds_table = preprocess_ann(arg.ann_path)
    end_time = timeit.default_timer()
    print('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()

    # read the variant file
    if mutation_mode == 'somatic_and_germline':
        mutation_dic_maf = parse_mutation_from_maf(maf_file_path)
        mutation_dic_vcf = parse_mutation_from_vcf(vcf_file_path)
    elif mutation_mode == 'germline':
        mutation_dic_maf = {} # empty dic
        mutation_dic_vcf = parse_mutation_from_vcf(vcf_file_path)
    elif mutation_mode == 'somatic':
        mutation_dic_maf = parse_mutation_from_maf(maf_file_path)
        mutation_dic_vcf = {} # empty dic
    elif mutation_mode == 'ref':
        mutation_dic_maf = {}
        mutation_dic_vcf = {}
    else:
        print('Mutation mode "%s" not recognized' % mutation_mode)
        sys.exit(1)

    # Personalized analysis
    # load splicegraph
    print('Loading splice graph ...')
    start_time = timeit.default_timer()
    with open(arg.splice_path, 'r') as graph_fp:
        (graph_data, graph_meta) = cPickle.load(graph_fp)  # both graph data and meta data
    end_time = timeit.default_timer()
    print('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()

    if arg.process_num == 0: # Default process all genes
        num = len(graph_data)
    else:
        num = arg.process_num

    # load graph metadata
    start_time = timeit.default_timer()
    if arg.count_path is not None:
        print('Loading count data ...')
        h5f = h5py.File(arg.count_path, 'r')
        seg_lookup_table, edge_lookup_table, strain_idx_table, segment_expr_info, edge_expr_info = parse_gene_metadata_info(h5f, arg.samples)
        end_time = timeit.default_timer()
        print('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
        print_memory_diags()
        # get the fp
        #strains = sp.array([x.split('.')[0] for x in h5f['strains'][:]])
        #size_factor = get_size_factor(strains, arg.libsize_path)
        size_factor = None
    else:
        seg_lookup_table = None
        edge_lookup_table = None
        strain_idx_table = None
        segment_expr_info = None
        edge_expr_info = None
        size_factor = None

    # read the intro of interest file gtex_junctions.hdf5
    junction_dict = parse_junction_meta_info(arg.gtex_junction_path)

    # process the genes according to the annotation file
    # add CDS starts and reading frames to the respective nodes
    print('Processing gene set ...')
    start_time = timeit.default_timer()
    anno_pickle = os.path.join(arg.output_dir, 'annotation_preprop.pickle')
    if os.path.exists(anno_pickle):
        print('...loading from preprocessed dump: %s' % anno_pickle)
        (graph_data, gene_cds_begin_dict) = cPickle.load(open(anno_pickle, 'r'))
    else:
        print('...computing from annotation')
        genes_preprocess(graph_data, gene_cds_begin_dict)
        cPickle.dump((graph_data, gene_cds_begin_dict), open(anno_pickle, 'w'), -1)
    end_time = timeit.default_timer()
    print('\tTime spent: {:.3f} seconds'.format(end_time - start_time))

    # process graph for each input sample
    for sample in arg.samples:
        # prepare for the output file
        output_path = os.path.join(arg.output_dir, sample)
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        peptide_file_path = os.path.join(output_path, mutation_mode + '_peptides.fa')
        meta_peptide_file_path = os.path.join(output_path, mutation_mode + '_metadata.tsv.gz')
        peptide_fp = open(peptide_file_path, 'w')
        meta_peptide_fp = gzip.open(meta_peptide_file_path, 'w')
        meta_header_line = "\t".join(['output_id','read_frame','gene_name', 'gene_chr', 'gene_strand','mutation_mode','peptide_weight','peptide_annotated',
                                    'junction_annotated','has_stop_codon','is_in_junction_list','is_isolated','variant_comb','seg_expr', 'exons_coor', 'vertex_idx','junction_expr'])
        meta_peptide_fp.write(meta_header_line + '\n')

        # go over each gene in splicegraph
        for gene_idx, gene in enumerate(graph_data[:num]):
            start_time = timeit.default_timer()
            print('%s %i/%i\n'%(sample, gene_idx, num))
            gene = graph_data[gene_idx]

            # Genes not contained in the annotation...
            if gene.name not in gene_cds_begin_dict or gene.name not in gene_to_transcript_table:
                gene.processed = False
                continue

            chrm = gene.chr.strip()
            if (sample, chrm) in mutation_dic_vcf.keys():
                mutation_sub_dict_vcf = mutation_dic_vcf[(sample, chrm)]
            else:
                mutation_sub_dict_vcf = None
            if (sample, chrm) in mutation_dic_maf.keys():
                mutation_sub_dict_maf = mutation_dic_maf[(sample, chrm)]
            else:
                mutation_sub_dict_maf = None

            if segment_expr_info is not None:
                sub_segment_expr_info = segment_expr_info[:, strain_idx_table[sample]]
                sub_edge_expr_info = edge_expr_info[:, strain_idx_table[sample]]
            else:
                sub_segment_expr_info = None
                sub_edge_expr_info = None

            if not junction_dict is None and chrm in junction_dict:
                junction_list = junction_dict[chrm]
            else:
                junction_list = None

            output_peptide_list, output_metadata_list = annotate_gene_opt(gene=gene,
                              ref_seq=seq_dict[chrm],
                              gene_idx=gene_idx,
                              seg_lookup_table=seg_lookup_table,
                              edge_lookup_table=edge_lookup_table,
                              size_factor=size_factor,
                              junction_list=junction_list,
                              segment_expr_info=sub_segment_expr_info,
                              edge_expr_info=sub_edge_expr_info,
                              transcript_to_cds_table=transcript_to_cds_table,
                              gene_to_transcript_table=gene_to_transcript_table,
                              mutation_mode=mutation_mode,
                              mutation_sub_dic_vcf=mutation_sub_dict_vcf,
                              mutation_sub_dic_maf=mutation_sub_dict_maf,
                              debug=arg.debug
                            )
            if arg.is_filter:
                output_metadata_list, output_peptide_list = get_filtered_output_list(output_metadata_list,output_peptide_list)
            meta_peptide_fp.write('\n'.join(output_metadata_list)+'\n')
            peptide_fp.write('\n'.join(output_peptide_list)+'\n')
            end_time = timeit.default_timer()
            print(gene_idx, end_time - start_time,'\n')


if __name__ == "__main__":
    main()
