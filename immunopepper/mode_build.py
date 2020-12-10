# Python libraries
""""""
# Core operation of ImmunoPepper. Traverse splicegraph and get kmer/peptide output
from collections import defaultdict
import h5py
import logging
import numpy as np
import os
import pathlib
import pickle
import signal as sig
import sys
import timeit


# immuno module
from .config import ExceptionWrapper
from .config import MyPool
from .io_ import collect_results
from .io_ import initialize_fp
from .io_ import remove_folder_list
from .io_ import save_backgrd_kmer_set
from .io_ import save_backgrd_pep_dict
from .io_ import save_gene_expr_distr
from .io_ import save_forgrd_kmer_dict
from .io_ import save_forgrd_pep_dict
from .mutations import get_mutation_mode_from_parser
from .mutations import get_sub_mutation_tuple
from .preprocess import genes_preprocess_all
from .preprocess import parse_junction_meta_info
from .preprocess import parse_gene_metadata_info
from .preprocess import preprocess_ann
from .traversal import collect_background_transcripts
from .traversal import collect_vertex_pairs
from .traversal import get_and_write_background_peptide_and_kmer
from .traversal import get_and_write_peptide_and_kmer
from .utils import check_chr_consistence
from .utils import create_libsize
from .utils import get_idx
from .utils import get_total_gene_expr
from .utils import print_memory_diags

### intermediate fix to load pickle files stored under previous version
from spladder.classes import gene as cgene
from spladder.classes import splicegraph as csplicegraph
from spladder.classes import segmentgraph as csegmentgraph
sys.modules['modules'] = cgene
sys.modules['modules.classes.gene'] = cgene
sys.modules['modules.classes.splicegraph'] = csplicegraph
sys.modules['modules.classes.segmentgraph'] = csegmentgraph
### end fix


def process_gene_batch_background(sample, genes, gene_idxs,  mutation , countinfo, genetable, arg, outbase, filepointer, compression=None, verbose=False):
    try:

        set_kmer_back =  defaultdict(set, {})
        dict_pept_backgrd = {}
        time_per_gene = []
        mem_per_gene = []
        all_gene_idxs = []

        for i, gene in enumerate(genes):
            ### measure time
            start_time = timeit.default_timer()

            # Genes not contained in the annotation...
            if gene.name not in genetable.gene_to_cds_begin or \
                    gene.name not in genetable.gene_to_ts:
                #logger.warning('>Gene name {} is not in the genetable and not processed, please check the annotation file.'.format(gene.name))
                continue

            idx = get_idx(countinfo, sample, gene_idxs[i])
            all_gene_idxs.append(gene_idxs[i])

            chrm = gene.chr.strip()
            sub_mutation = get_sub_mutation_tuple(mutation, sample, chrm)
            ref_mut_seq = collect_background_transcripts(gene=gene, ref_seq_file=arg.ref_path, chrm=chrm, mutation=sub_mutation)

            # Gene counts information
            if countinfo:
                gidx = countinfo.gene_idx_dict[gene.name]
                count_segments = np.arange(countinfo.gene_id_to_segrange[gidx][0], countinfo.gene_id_to_segrange[gidx][1])
                with h5py.File(countinfo.h5fname, 'r') as h5f:
                    seg_counts = h5f['segments'][count_segments, idx.sample]

            get_and_write_background_peptide_and_kmer(peptide_dict = dict_pept_backgrd,
                                                      kmer_dict = set_kmer_back,
                                                      gene=gene,
                                                      ref_mut_seq=ref_mut_seq,
                                                       gene_table=genetable,
                                                       countinfo=countinfo,
                                                       seg_counts=seg_counts,
                                                       Idx=idx,
                                                       kmer=arg.kmer)

            time_per_gene.append(timeit.default_timer() - start_time)
            mem_per_gene.append(print_memory_diags(disable_print=True))

        save_backgrd_pep_dict(dict_pept_backgrd, filepointer, compression, outbase, verbose)
        dict_pept_backgrd.clear()
        for kmer_length in set_kmer_back:
            save_backgrd_kmer_set(set_kmer_back[kmer_length], filepointer, kmer_length, compression, outbase, verbose)
        set_kmer_back.clear()

        if all_gene_idxs:
            logging.info("> {}: annotation graph {}/{} processed, max time cost: {}, memory cost:{} GB for gene batch".format(sample,
                                                                                      all_gene_idxs[-1]  + 1,
                                                                                      len(gene_id_list),
                                                                                      np.max(time_per_gene),
                                                                                      np.max(mem_per_gene)))
        exception_ = None

    except Exception as e:
        exception_ = ExceptionWrapper(e)

    return exception_



def process_gene_batch_foreground(sample, graph_samples, genes, genes_info, gene_idxs, total_genes, all_read_frames, mutation, junction_dict, countinfo, genetable, arg, outbase, filepointer, compression, verbose):
    try:
        pathlib.Path(outbase).mkdir(exist_ok=True, parents=True)
        dict_kmer_foregr = defaultdict(dict, {})
        dict_pept_forgrd = {}
        time_per_gene = []
        mem_per_gene = []
        all_gene_idxs = []
        gene_expr = []
        complexity_cap = 20000

        test = ['sample', 'graph_samples','genes', 'genes_info', 'gene_idxs', 'total_genes', 'all_read_frames', 'mutation', 'junction_dict',
         'countinfo', 'genetable', 'arg', 'outbase', 'filepointer', 'compression', 'verbose']
        for i, gene in enumerate(genes):
            ### measure time
            start_time = timeit.default_timer()

            # Genes not contained in the annotation...
            if gene.name not in genetable.gene_to_cds_begin or \
               gene.name not in genetable.gene_to_ts:
                #logging.warning('>Gene {} is not in the genetable and not processed, please check the annotation file.'.format(gene.name))
                continue

            # Genes with highly complex splicegraphs
            if (len(gene.splicegraph.vertices[1]) > complexity_cap):
                logging.warning('>Gene {} has a edge complexity > {}, not processed'.format(gene.name, complexity_cap))
                continue

            idx = get_idx(countinfo, sample, gene_idxs[i])
            logging.info("gene processed {}".format(i))
            # Gene counts information
            if countinfo:
                gidx = countinfo.gene_idx_dict[gene.name]

                with h5py.File(countinfo.h5fname, 'r') as h5f:
                        if countinfo.gene_idx_dict[gene.name] not in countinfo.gene_id_to_edgerange or \
                                countinfo.gene_idx_dict[gene.name] not in countinfo.gene_id_to_segrange:
                            edge_idxs = None
                            edge_counts = None

                        else:
                            edge_gene_idxs = np.arange(countinfo.gene_id_to_edgerange[gidx][0], countinfo.gene_id_to_edgerange[gidx][1])
                            edge_idxs = h5f['edge_idx'][list(edge_gene_idxs)].astype('int')
                            if arg.cross_graph_expr:
                                edge_counts = h5f['edges'][edge_gene_idxs,:] # will compute expression on whole graph
                            else:
                                edge_counts = h5f['edges'][edge_gene_idxs,idx.sample]
                        seg_gene_idxs = np.arange(countinfo.gene_id_to_segrange[gidx][0],
                                                  countinfo.gene_id_to_segrange[gidx][1])
                        if arg.cross_graph_expr:
                            seg_counts = h5f['segments'][seg_gene_idxs, :]
                        else:
                            seg_counts = h5f['segments'][seg_gene_idxs, idx.sample]

            if not arg.cross_graph_expr: #TODO deal with gene_expression in sample in this case
                gene_expr.append((gene.name, get_total_gene_expr(gene, countinfo, idx, seg_counts)))

            chrm = gene.chr.strip()
            sub_mutation = get_sub_mutation_tuple(mutation, sample, chrm)
            junction_list = None
            if not junction_dict is None and chrm in junction_dict:
                junction_list = junction_dict[chrm]

            vertex_pairs, \
            ref_mut_seq, \
            exon_som_dict = collect_vertex_pairs(gene=gene,
                                                 gene_info=genes_info[i],
                                                 ref_seq_file=arg.ref_path, #seq_dict[chrm],
                                                 chrm=chrm,
                                                 idx=idx,
                                                 mutation=sub_mutation,
                                                 all_read_frames=all_read_frames,
                                                 disable_concat=arg.disable_concat,
                                                 kmer=arg.kmer,
                                                 filter_redundant=arg.filter_redundant)

            get_and_write_peptide_and_kmer(peptide_dict = dict_pept_forgrd,
                                            kmer_dict = dict_kmer_foregr,
                                            gene=gene,
                                            all_vertex_pairs=vertex_pairs,
                                            ref_mut_seq=ref_mut_seq,
                                            idx=idx,
                                            exon_som_dict=exon_som_dict,
                                            countinfo=countinfo,
                                            mutation=sub_mutation,
                                            table=genetable,
                                            size_factor=None,
                                            junction_list=junction_list,
                                            kmer=arg.kmer,
                                            output_silence=arg.output_silence,
                                            outbase=outbase,
                                            edge_idxs=edge_idxs,
                                            edge_counts=edge_counts,
                                            seg_counts=seg_counts,
                                            cross_graph_expr=arg.cross_graph_expr,
                                            filepointer=filepointer,
                                            graph_samples=graph_samples,
                                            verbose_save=verbose
            )


            time_per_gene.append(timeit.default_timer() - start_time)
            mem_per_gene.append(print_memory_diags(disable_print=True))
            all_gene_idxs.append(gene_idxs[i])
        
        save_gene_expr_distr(gene_expr, filepointer, outbase, compression, verbose)
        save_forgrd_pep_dict(dict_pept_forgrd, filepointer, compression, outbase, arg.output_fasta, verbose)
        dict_pept_forgrd.clear()
        if not arg.cross_graph_expr:
            for kmer_length in dict_kmer_foregr:
                save_forgrd_kmer_dict(dict_kmer_foregr[kmer_length], filepointer, kmer_length, compression, outbase, verbose)
            dict_kmer_foregr.clear()
        if arg.cross_graph_expr and filepointer.kmer_segm_expr_fp['pqwriter'] is not None:
            filepointer.kmer_segm_expr_fp['pqwriter'].close()
            filepointer.kmer_edge_expr_fp['pqwriter'].close()

        if all_gene_idxs:
            logging.info("> {}: sample graph {}/{} processed, max time cost: {}, memory cost:{} GB for gene batch".format(sample,
                                                                                          all_gene_idxs[-1]  + 1,
                                                                                          len(gene_id_list),
                                                                                          np.max(time_per_gene),
                                                                                        np.max(mem_per_gene)))
        exception_ = None

    
    except Exception as e:

        exception_ = ExceptionWrapper(e)
    return exception_




def mode_build(arg):
    # read and process the annotation file
    logging.info(">>>>>>>>> Build: Start Preprocessing")
    logging.info('Building lookup structure ...')
    start_time = timeit.default_timer()
    genetable,chromosome_set = preprocess_ann(arg.ann_path)
    end_time = timeit.default_timer()
    logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()

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
    #graph_data = graph_data[[3170]] #TODO remove
    #graph_data = graph_data[400:5400]
    #graph_data = graph_data[0:20]
    all_read_frames = arg.all_read_frames

    check_chr_consistence(chromosome_set,mutation,graph_data)

    if arg.process_num == 0:  # Default process all genes
        num = len(graph_data)
    else:
        num = arg.process_num
    
    # load graph metadata
    start_time = timeit.default_timer()
    if arg.count_path is not None:
        logging.info('Loading count data ...')
        countinfo, graph_samples = parse_gene_metadata_info(arg.count_path, arg.samples, arg.cross_graph_expr)
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
    graph_info = genes_preprocess_all(graph_data, genetable.gene_to_cds_begin, arg.parallel, all_read_frames)
    end_time = timeit.default_timer()
    logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()
    logging.info(">>>>>>>>> Finish Preprocessing")
    expr_distr_dict = {}
    # process graph for each input sample
    output_libszie_fp = os.path.join(arg.output_dir,'expression_counts.libsize.tsv')
    logging.info(">>>>>>>>> Start traversing splicegraph")

    #if arg.parallel > 1:
    global gene_id_list
    global sample
    global filepointer


    for sample in arg.samples:
        logging.info(">>>> Processing sample {}, there are {} graphs in total".format(sample,num))

        # prepare for the output file
        output_path = os.path.join(arg.output_dir, sample)
        if not os.path.isdir(output_path):
            os.makedirs(output_path)


        if arg.compressed:
            gzip_tag = '.gz'
            pq_compression = 'GZIP'
        else:
            pq_compression = None
            gzip_tag = ''

        filepointer = initialize_fp(output_path, mutation.mode, gzip_tag,
                  arg.kmer, arg.output_fasta, arg.cross_graph_expr)

        # go over each gene in splicegraph
        gene_id_list = list(range(0,num))

        if arg.parallel > 1:

            logging.info('Parallel: {} Threads'.format(arg.parallel))
            batch_size = min(num, arg.batch_size)
            verbose_save = False
#            # Build the background
#            logging.info(">>>>>>>>> Start Background processing")
#            pool_f = MyPool(processes=arg.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
#            for i in range(0, len(gene_id_list), batch_size):
#                gene_idx = gene_id_list[i:min(i + batch_size, len(gene_id_list))]
#                outbase = os.path.join(output_path, 'tmp_out_{}_{}'.format(arg.mutation_mode, i))
#                pathlib.Path(outbase).mkdir(exist_ok= True, parents= True)
#
#                _ = pool_f.submit(process_gene_batch_background, (sample, graph_data[gene_idx], gene_idx, mutation, countinfo, genetable, arg, outbase, filepointer, None, verbose_save))
#            pool_f.terminate()
#
            # Build the foreground and remove the background if needed
            logging.info(">>>>>>>>> Start Foreground processing")
            pool_f = MyPool(processes=arg.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))

            gene_batches = [(i, gene_id_list[i:min(i + batch_size, len(gene_id_list))]) for i in range(0, len(gene_id_list), batch_size)]

            args = [(sample, graph_samples, graph_data[gene_idx], graph_info[gene_idx], gene_idx, len(
                gene_id_list), all_read_frames, mutation, junction_dict, countinfo, genetable, arg,
                  os.path.join(output_path, 'tmp_out_{}_{}'.format(arg.mutation_mode, i)), filepointer, None, verbose_save) for i, gene_idx in gene_batches ]

            result = pool_f.submit(process_gene_batch_foreground, args)
            pool_f.terminate()



            # Collects and pools the files of each batch
            logging.debug('start collecting results')
            collect_results(filepointer.background_peptide_fp, output_path, pq_compression, arg.mutation_mode)
            collect_results(filepointer.junction_peptide_fp, output_path, pq_compression, arg.mutation_mode)
            collect_results(filepointer.junction_meta_fp, output_path, pq_compression, arg.mutation_mode)
            collect_results(filepointer.junction_kmer_fp, output_path, pq_compression, arg.mutation_mode, arg.kmer)
            collect_results(filepointer.background_kmer_fp, output_path, pq_compression, arg.mutation_mode, arg.kmer)
            collect_results(filepointer.gene_expr_fp, output_path, pq_compression, arg.mutation_mode)
            collect_results(filepointer.kmer_segm_expr_fp, output_path, pq_compression, arg.mutation_mode)
            collect_results(filepointer.kmer_edge_expr_fp, output_path, pq_compression, arg.mutation_mode)
#            remove_folder_list(os.path.join(output_path, 'tmp_out_{}'.format(arg.mutation_mode)))

        else:
            logging.info('Not Parallel')
            # Build the background
            logging.info(">>>>>>>>> Start Background processing")
            process_gene_batch_background(sample, graph_data, gene_id_list, mutation, countinfo, genetable, arg, output_path, filepointer, pq_compression, verbose=True)
            # Build the foreground and remove the background if needed
            logging.info(">>>>>>>>> Start Foreground processing")
            process_gene_batch_foreground( sample, graph_samples, graph_data, graph_info, gene_id_list, len(gene_id_list), all_read_frames, mutation, junction_dict,
                             countinfo, genetable, arg, output_path, filepointer, pq_compression, verbose=True)

        if not arg.cross_graph_expr:
            create_libsize(filepointer.gene_expr_fp,output_libszie_fp, sample)




