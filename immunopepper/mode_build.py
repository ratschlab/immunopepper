# Python libraries
""""""
# Core operation of ImmunoPepper. Traverse splicegraph and get kmer/peptide output
import logging
import gzip
import os
import pickle
import sys
import timeit
import multiprocessing as mp
import signal as sig
from functools import partial

# External libraries
import Bio.SeqIO as BioIO
import gc
import h5py
import numpy as np

# immuno module
from .io import convert_namedtuple_to_str
from .io import gz_and_normal_open
from .io import load_pickled_graph
from .io import print_memory_diags
from .io import write_gene_expr
from .io import write_namedtuple_list
from .mutations import get_mutation_mode_from_parser
from .mutations import get_sub_mutation_tuple
from .namedtuples import Filepointer
from .preprocess import genes_preprocess_all
from .preprocess import parse_junction_meta_info
from .preprocess import parse_gene_metadata_info
from .preprocess import preprocess_ann
from .traversal import collect_vertex_pairs
from .traversal import get_and_write_background_peptide_and_kmer
from .traversal import get_and_write_peptide_and_kmer
from .utils import check_chr_consistence
from .utils import create_libsize
from .utils import get_idx
from .utils import get_total_gene_expr

### intermediate fix to load pickle files stored under previous version
from spladder.classes import gene as cgene
from spladder.classes import splicegraph as csplicegraph
from spladder.classes import segmentgraph as csegmentgraph
sys.modules['modules'] = cgene
sys.modules['modules.classes.gene'] = cgene
sys.modules['modules.classes.splicegraph'] = csplicegraph
sys.modules['modules.classes.segmentgraph'] = csegmentgraph
### end fix


def process_single_gene(sample, gene, gene_info, gene_idx, mutation, junction_dict, countinfo, genetable, arg):
    ### set result dict
    R = dict()
    R['gene_name'] = gene.name

    # Genes not contained in the annotation...
    if gene.name not in genetable.gene_to_cds_begin or gene.name not in genetable.gene_to_ts:
        gene.processed = False
        logging.warning('>Gene name {} is not in the genetable and not processed, please check the annotation file.'.format(gene.name))
        R['processed'] = False
        return R
    R['processed'] = True

    idx = get_idx(countinfo, sample, gene_idx)
    R['total_expr'] = get_total_gene_expr(gene, countinfo, idx)

    chrm = gene.chr.strip()
    sub_mutation = get_sub_mutation_tuple(mutation, sample, chrm)
    junction_list = None
    if not junction_dict is None and chrm in junction_dict:
        junction_list = junction_dict[chrm]

    vertex_pairs, ref_mut_seq, exon_som_dict = collect_vertex_pairs(gene=gene,
                                                                    gene_info=gene_info,
                                                                    ref_seq_file=arg.ref_path, #seq_dict[chrm],
                                                                    chrm=chrm,
                                                                    idx=idx,
                                                                    mutation=sub_mutation,
                                                                    disable_concat=arg.disable_concat,
                                                                    kmer=arg.kmer,
                                                                    filter_redundant=arg.filter_redundant)

    R['background_peptide_list'], R['background_kmer_lists'] = get_and_write_background_peptide_and_kmer(gene=gene, 
                                                                    ref_mut_seq=ref_mut_seq,
                                                                    gene_table=genetable,
                                                                    countinfo=countinfo, 
                                                                    Idx=idx,
                                                                    kmer=arg.kmer)
    R['output_metadata_list'], R['output_peptide_list'], R['output_kmer_lists'] = get_and_write_peptide_and_kmer(gene=gene, 
                                                                               vertex_pairs=vertex_pairs, 
                                                                               background_pep_list=R['background_peptide_list'],
                                                                               ref_mut_seq=ref_mut_seq, 
                                                                               idx=idx, 
                                                                               exon_som_dict=exon_som_dict,
                                                                               countinfo=countinfo,
                                                                               mutation=sub_mutation,
                                                                               table=genetable,
                                                                               size_factor=None,
                                                                               junction_list=junction_list, 
                                                                               kmer=arg.kmer,
                                                                               output_silence=arg.output_silence)

    return R


def write_gene_result(gene_result, filepointer):

    ### define fields relevant for output
    back_pep_field_list = ['id', 'new_line', 'peptide']
    for peptide in gene_result['background_peptide_list']:
        filepointer.background_peptide_fp.write(convert_namedtuple_to_str(peptide, back_pep_field_list) + '\n')

    back_kmer_field_list = ['kmer', 'id', 'expr', 'is_cross_junction']
    for kmer_list in gene_result['background_kmer_lists']:
        write_namedtuple_list(filepointer.background_kmer_fp, kmer_list, back_kmer_field_list)

    junc_pep_field_list = ['output_id', 'id', 'new_line', 'peptide']
    for peptide in gene_result['output_peptide_list']:
        filepointer.junction_peptide_fp.write(convert_namedtuple_to_str(peptide, junc_pep_field_list) + '\n')

    meta_field_list = ['output_id', 'read_frame', 'gene_name', 'gene_chr', 'gene_strand', 'mutation_mode', 'peptide_annotated',
                       'junction_annotated', 'has_stop_codon', 'is_in_junction_list', 'is_isolated', 'variant_comb',
                       'variant_seg_expr', 'modified_exons_coord','original_exons_coord', 'vertex_idx', 'junction_expr', 'segment_expr']
    for output_metadata in gene_result['output_metadata_list']:
        filepointer.junction_meta_fp.write(convert_namedtuple_to_str(output_metadata, meta_field_list) + '\n')

    kmer_field_list = ['kmer', 'id', 'expr', 'is_cross_junction', 'junction_count']
    for kmer_list in gene_result['output_kmer_lists']:
        write_namedtuple_list(filepointer.junction_kmer_fp, kmer_list, kmer_field_list)


def mode_build(arg):
    # read and process the annotation file
    logging.info(">>>>>>>>> Build: Start Preprocessing")
    logging.info('Building lookup structure ...')
    start_time = timeit.default_timer()
    genetable,chromosome_set = preprocess_ann(arg.ann_path)
    end_time = timeit.default_timer()
    logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()

    # load genome sequence data
    #seq_dict = {}
    #start_time = timeit.default_timer()
    #logging.info('Parsing genome sequence ...')
    #for record in BioIO.parse(arg.ref_path, "fasta"):
    #    if record.id in chromosome_set:
    #        seq_dict[record.id] = str(record.seq).strip()
    #end_time = timeit.default_timer()
    #logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    #print_memory_diags()

    # read the variant file
    mutation = get_mutation_mode_from_parser(arg)
    # load splicegraph
    logging.info('Loading splice graph ...')
    start_time = timeit.default_timer()
    with open(arg.splice_path, 'rb') as graph_fp:
        (graph_data, graph_meta) = pickle.load(graph_fp, encoding='latin1')  # both graph data and meta data
    end_time = timeit.default_timer()
    logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()
    
    ### DEBUG
    #graph_data = graph_data[[3170]]
    #graph_data = graph_data[:1000]

    check_chr_consistence(chromosome_set,mutation,graph_data)

    if arg.process_num == 0:  # Default process all genes
        num = len(graph_data)
    else:
        num = arg.process_num

    # load graph metadata
    start_time = timeit.default_timer()
    if arg.count_path is not None:
        logging.info('Loading count data ...')
        countinfo = parse_gene_metadata_info(arg.count_path, arg.samples)
        end_time = timeit.default_timer()
        logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
        print_memory_diags()
        #size_factor = get_size_factor(strains, arg.libsize_path)
        #size_factor = None
    else:
        countinfo = None
        size_factor = None

    # read the intron of interest file gtex_junctions.hdf5
    junction_dict = parse_junction_meta_info(arg.gtex_junction_path)

    # process the genes according to the annotation file
    # add CDS starts and reading frames to the respective nodes
    logging.info('Add reading frame to splicegraph ...')
    start_time = timeit.default_timer()
    graph_info = genes_preprocess_all(graph_data, genetable.gene_to_cds_begin, arg.parallel)
    end_time = timeit.default_timer()
    logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()
    logging.info(">>>>>>>>> Finish Preprocessing")
    expr_distr_dict = {}
    # process graph for each input sample
    output_libszie_fp = os.path.join(arg.output_dir,'expression_counts.libsize.tsv')
    logging.info(">>>>>>>>> Start traversing splicegraph")

    if arg.parallel > 1:
        global cnt
        global gene_name_expr_distr
        global expr_distr
        global filepointer
        global gene_id_list

    for sample in arg.samples:
        logging.info(">>>> Processing sample {}, there are {} graphs in total".format(sample,num))
        error_gene_num = 0
        time_list = []
        memory_list = []
        expr_distr = []
        gene_name_expr_distr = []
        # prepare for the output file
        output_path = os.path.join(arg.output_dir, sample)
        if not os.path.isdir(output_path):
            os.makedirs(output_path)

        junction_meta_file_path = os.path.join(output_path, mutation.mode + '_metadata.tsv.gz')
        meta_peptide_fp = gz_and_normal_open(junction_meta_file_path,'w')
        gzip_tag = ''
        if arg.compressed:
            gzip_tag = '.gz'
        junction_peptide_file_path = os.path.join(output_path, mutation.mode + '_peptides.fa'+gzip_tag)
        background_peptide_file_path = os.path.join(output_path, mutation.mode + '_back_peptides.fa'+gzip_tag)
        junction_kmer_file_path = os.path.join(output_path, mutation.mode + '_junction_kmer.txt'+gzip_tag)
        background_kmer_file_path = os.path.join(output_path, mutation.mode + '_back_kmer.txt'+gzip_tag)
        gene_expr_file_path = os.path.join(output_path, 'gene_expression_detail.tsv'+gzip_tag)
        peptide_fp = gz_and_normal_open(junction_peptide_file_path, 'w')
        background_fp = gz_and_normal_open(background_peptide_file_path, 'w')
        junction_kmer_fp = gz_and_normal_open(junction_kmer_file_path, 'w')
        background_kmer_fp = gz_and_normal_open(background_kmer_file_path, 'w')
        gene_expr_fp = gz_and_normal_open(gene_expr_file_path, 'w')

        filepointer = Filepointer(peptide_fp, meta_peptide_fp, background_fp, junction_kmer_fp, background_kmer_fp)

        meta_field_list = ['output_id', 'read_frame', 'gene_name', 'gene_chr', 'gene_strand', 'mutation_mode',
                           'peptide_annotated',
                           'junction_annotated', 'has_stop_codon', 'is_in_junction_list', 'is_isolated', 'variant_comb',
                           'variant_seg_expr',
                           'modified_exons_coord', 'original_exons_coord', 'vertex_idx', 'junction_expr',
                           'segment_expr']
        meta_peptide_fp.write(('\t'.join(meta_field_list) + '\n'))

        junction_kmer_field_list = ['kmer','gene_name','seg_expr','is_crossjunction','junction_expr']
        junction_kmer_fp.write(('\t'.join(junction_kmer_field_list)+'\n'))

        background_kmer_field_list =  ['kmer','gene_name','seg_expr','is_crossjunction']
        background_kmer_fp.write(('\t'.join(background_kmer_field_list)+'\n'))
        # go over each gene in splicegraph
        gene_id_list = list(range(0,num))
        if arg.parallel > 1:
            cnt = 0

            def process_result(gene_result):
                global cnt
                global gene_name_expr_distr
                global expr_distr
                global filepointer
                global gene_id_list
                logging.info('> Finished processing Gene {} ({}/{})'.format(gene_result['gene_name'], cnt, len(gene_id_list)))
                cnt += 1
                if gene_result['processed']:
                    gene_name_expr_distr.append((gene_result['gene_name'], gene_result['total_expr']))
                    expr_distr.append(gene_result['total_expr'])
                    write_gene_result(gene_result, filepointer)
                del gene_result

            pool = mp.Pool(processes=arg.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
            for gene_idx in gene_id_list:
                _ = pool.apply_async(process_single_gene, args=(sample, graph_data[gene_idx], graph_info[gene_idx], gene_idx, mutation, junction_dict, countinfo, genetable, arg,), callback=process_result)
            pool.close()
            pool.join()
        else:
            for gene_idx in gene_id_list:
                gene = graph_data[gene_idx]
                gene_info = graph_info[gene_idx]
                start_time = timeit.default_timer()
                gene_result = process_single_gene(sample, gene, gene_info, gene_idx, mutation, junction_dict, countinfo, genetable, arg)
                if gene_result['processed']:
                    gene_name_expr_distr.append((gene.name, gene_result['total_expr']))
                    expr_distr.append(gene_result['total_expr'])
                    write_gene_result(gene_result, filepointer)
                end_time = timeit.default_timer()
                memory = print_memory_diags(disable_print=True)
                logging.info(">{}: {}/{} processed, time cost: {}, memory cost:{} GB ".format(sample,gene_idx+1, len(gene_id_list),end_time - start_time,memory))
                time_list.append(end_time - start_time)
                memory_list.append(memory)

                if len(arg.samples) == 1:
                    graph_data[gene_idx] = None
                    gc.collect()
        if memory_list and time_list:
            max_memory,max_time = max(memory_list),max(time_list)
            max_memory_id,max_time_id = np.argmax(memory_list),np.argmax(time_list)
            if num-error_gene_num > 0:
                mean_memory, mean_time = sum(memory_list)/(num-error_gene_num),sum(time_list)/(num-error_gene_num)
            else:
                mean_memory, mean_time = 0,0
            logging.info(">>>> Finish sample {}. Errors existed in {}/{} genes. Might Need further check. "
                         "Max memory :{} GB, Max time :{} seconds, Max memory gene ID:{}, Max time gene ID:{}, "
                         "Average memory cost:{} GB, Average time cost:{} seconds".format(sample,error_gene_num,num,max_memory,max_time,max_memory_id,max_time_id,
                                                                               mean_memory,mean_time))
        else:
            logging.info(">>>> No gene is processed during this running.")
        expr_distr_dict[sample] = expr_distr
        write_gene_expr(gene_expr_fp,gene_name_expr_distr)
        create_libsize(expr_distr_dict,output_libszie_fp)
    
        ### close files
        filepointer.junction_peptide_fp.flush()
        filepointer.junction_peptide_fp.close()
        filepointer.junction_meta_fp.flush()
        filepointer.junction_meta_fp.close()
        filepointer.background_peptide_fp.flush()
        filepointer.background_peptide_fp.close()
        filepointer.junction_kmer_fp.flush()
        filepointer.junction_kmer_fp.close()
        filepointer.background_kmer_fp.flush()
        filepointer.background_kmer_fp.close()
    logging.info(">>>>>>>>> Build: Finish traversing splicegraph in mutation mode {}.\n".format(mutation.mode))

