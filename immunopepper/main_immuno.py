

# Python libraries
import sys
import pickle
import os
import timeit
import argparse
import gzip
import logging

# External libraries
import Bio.SeqIO as BioIO
import h5py

# immuno module
from immunopepper.immuno_print import print_memory_diags
from immunopepper.immuno_preprocess import genes_preprocess,preprocess_ann,parse_gene_metadata_info,parse_junction_meta_info
from immunopepper.immuno_mutation import get_mutation_mode_from_parser,get_sub_mutation_tuple
from immunopepper.immuno_model import get_simple_metadata, get_and_write_peptide_and_kmer,get_and_write_background_peptide_and_kmer
from immunopepper.immuno_nametuple import Option, Filepointer
from immunopepper.io_utils import load_pickled_graph
from immunopepper.utils import get_idx,create_libsize

def parse_arguments(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", nargs='+', help="the sample names, can specify more than one sample", required=False, default='')
    parser.add_argument("--output_dir", help="specify the output directory [default: tests]", required=False, default='test')
    parser.add_argument("--ann_path", help="specify the absolute path of annotation file", required=False)
    parser.add_argument("--splice_path", help="specify the absolute path of splicegraph file", required=False)
    parser.add_argument("--ref_path", help="specify the absolute path of reference gene file to the work_dir", required=False)
    parser.add_argument("--libsize_path", nargs='?', help="specify the absolute path to expression library sizes",required=False, default=None)
    parser.add_argument("--vcf_path", help="specify the absolute path of vcf file", required=False, default='')
    parser.add_argument("--maf_path", help="specify the absolute path of maf file", required=False, default='')
    parser.add_argument("--count_path",help="specify the absolute path of the count h5 file", required=False, default=None)
    parser.add_argument("--gtex_junction_path",help="specify the absolute path the the gtex_junction h5 file", required=False, default=None)
    parser.add_argument("--process_num", type=int, help="Only process the first *process_num* gene in the splicegraph,default,0, means process all", required=False, default=0)
    parser.add_argument("--filter_redundant", help="apply redundancy filter to the exon list", action="store_true", required=False, default=False)
    parser.add_argument("--debug", help="generate debug output", action="store_true", required=False, default=False)
    parser.add_argument("--mutation_mode", help="specify the mutation mdoe", required=False, default='ref')
    parser.add_argument("--heter_code", type=int, help="if count expression data is provided in h5 format, specify the code for heterzygous", default=0)
    parser.add_argument("--kmer", type=int, help="specify the k for kmer output", required=False, default=0)
    parser.add_argument("--output_silence",help="output mutated peptide even it is the same as reference peptide", action="store_true",default=False)
    if len(argv) < 2:
        parser.print_help()
        sys.exit(1)

    pargs = parser.parse_args(argv)
    return pargs

def main(arg):
    # load genome sequence data
    seq_dict = {}
    start_time = timeit.default_timer()
    interesting_chr = list(map(str, list(range(1, 23)))) + ["X", "Y", "MT"]
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
    genetable = preprocess_ann(arg.ann_path)
    end_time = timeit.default_timer()
    print('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()

    # read the variant file
    mutation = get_mutation_mode_from_parser(arg)

    # load splicegraph
    print('Loading splice graph ...')
    start_time = timeit.default_timer()
    with open(arg.splice_path, 'rb') as graph_fp:
        (graph_data, graph_meta) = pickle.load(graph_fp, encoding='latin1')  # both graph data and meta data
    end_time = timeit.default_timer()
    print('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()

    if arg.process_num == 0:  # Default process all genes
        num = len(graph_data)
    else:
        num = arg.process_num

    # load graph metadata
    start_time = timeit.default_timer()
    if arg.count_path is not None:
        print('Loading count data ...')
        h5f = h5py.File(arg.count_path, 'r')
        countinfo = parse_gene_metadata_info(h5f, arg.samples)
        edges, segments, sample_idx_table = countinfo.edges,countinfo.segments,countinfo.sample_idx_table
        end_time = timeit.default_timer()
        print('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
        print_memory_diags()
        # get the fp
        #strains = sp.array([x.split('.')[0] for x in h5f['strains'][:]])
        #size_factor = get_size_factor(strains, arg.libsize_path)
        #size_factor = None
    else:
        segments = None
        edges = None
        sample_idx_table = None
        size_factor = None

    # read the intron of interest file gtex_junctions.hdf5
    junction_dict = parse_junction_meta_info(arg.gtex_junction_path)

    # process the genes according to the annotation file
    # add CDS starts and reading frames to the respective nodes
    print('Processing gene set ...')
    start_time = timeit.default_timer()
    print('...computing from annotation')
    genes_preprocess(graph_data, genetable.gene_to_cds_begin)
    end_time = timeit.default_timer()
    print('\tTime spent: {:.3f} seconds'.format(end_time - start_time))

    expr_distr_dict = {}
    # process graph for each input sample
    output_libszie_fp = os.path.join(arg.output_dir,'expression_counts.libsize.tsv')
    option = Option(output_silence=arg.output_silence,debug=arg.debug,filter_redundant=arg.filter_redundant,kmer=arg.kmer)
    for sample in arg.samples:
        expr_distr = []
        # prepare for the output file
        output_path = os.path.join(arg.output_dir, sample)
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        log_dir = os.path.join(output_path, 'error.log')
        logging.basicConfig(level=logging.DEBUG,
                            filename=log_dir, filemode="a+",
                            format="%(asctime)-15s %(levelname)-8s %(message)s")

        junction_peptide_file_path = os.path.join(output_path, mutation.mode + '_peptides.fa')
        junction_meta_file_path = os.path.join(output_path, mutation.mode + '_metadata.tsv.gz')
        background_peptide_file_path = os.path.join(output_path, mutation.mode + '_back_peptides.fa')
        junction_kmer_file_path = os.path.join(output_path, mutation.mode + '_junction_kmer.txt')
        background_kmer_file_path = os.path.join(output_path, mutation.mode + '_back_kmer.txt')

        peptide_fp = open(junction_peptide_file_path, 'w')
        meta_peptide_fp = gzip.open(junction_meta_file_path, 'wt')
        background_fp = open(background_peptide_file_path,'w')
        junction_kmer_fp = open(junction_kmer_file_path, 'w')
        background_kmer_fp = open(background_kmer_file_path, 'w')

        filepointer = Filepointer(peptide_fp,meta_peptide_fp,background_fp,junction_kmer_fp,background_kmer_fp)

        meta_field_list = ['output_id','read_frame','gene_name', 'gene_chr', 'gene_strand','mutation_mode','peptide_weight','peptide_annotated',
                                    'junction_annotated','has_stop_codon','is_in_junction_list','is_isolated','variant_comb','variant_seg_expr',
                                      'exons_coor', 'vertex_idx','junction_expr','segment_expr']
        meta_peptide_fp.write(('\t'.join(meta_field_list) + '\n'))

        # go over each gene in splicegraph
        gene_id_list = list(range(0,num))
        for gene_idx in gene_id_list:
            #if 1:
            try:
                gene = graph_data[gene_idx]
                start_time = timeit.default_timer()
                print('%s %i/%i\n'%(sample, gene_idx, num))
                idx = get_idx(sample_idx_table,sample,gene_idx)
                # Genes not contained in the annotation...
                if gene.name not in genetable.gene_to_cds_begin or gene.name not in genetable.gene_to_ts:
                    gene.processed = False
                    continue

                chrm = gene.chr.strip()
                sub_mutation = get_sub_mutation_tuple(mutation,sample, chrm)
                if not junction_dict is None and chrm in junction_dict:
                    junction_list = junction_dict[chrm]
                else:
                    junction_list = None

                final_simple_meta, ref_mut_seq, exon_som_dict = get_simple_metadata(gene=gene,ref_seq=seq_dict[chrm],
                                                                                    idx=idx,mutation=sub_mutation,option=option)

                background_pep_list = get_and_write_background_peptide_and_kmer(gene=gene,ref_mut_seq=ref_mut_seq,table=genetable,
                                                          Idx=idx,filepointer=filepointer,option=option,Segments=segments)

                total_expr = get_and_write_peptide_and_kmer(gene=gene, final_simple_meta=final_simple_meta, background_pep_list=background_pep_list,
                                               ref_mut_seq=ref_mut_seq, idx=idx, exon_som_dict=exon_som_dict, segments=segments,
                                               edges=edges, mutation=sub_mutation,table=genetable,option=option,size_factor=None,
                                               junction_list=junction_list, filepointer=filepointer)
                expr_distr.append(total_expr)
                end_time = timeit.default_timer()
                print(gene_idx, end_time - start_time,'\n')

            except Exception as e:
                # should also print the error
                logging.exception("Exception occured in gene %d, %s mode, sample %s " % (gene_idx,arg.mutation_mode,sample))

        expr_distr_dict[sample] = expr_distr
    create_libsize(expr_distr_dict,output_libszie_fp)

def cmd_entry():
    arg = parse_arguments(sys.argv[1:])
    main(arg)


if __name__ == "__main__":
    cmd_entry()
