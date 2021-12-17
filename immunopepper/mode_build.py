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
from .preprocess import parse_gene_choices
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

def mapper_funct(tuple_arg):
    process_gene_batch_foreground(*tuple_arg)

def mapper_funct_back(tuple_arg):
    process_gene_batch_background(*tuple_arg)

def process_gene_batch_background(output_sample, mutation_sample, genes, gene_idxs, all_read_frames, mutation, countinfo_b, genetable, arg, outbase, filepointer, compression=None, verbose=False):
    if arg.parallel > 1:
        batch_name = int(outbase.split('/')[-1].split('_')[-1])
    else:
        batch_name = 'all'

    if (arg.parallel==1) or (not os.path.exists(os.path.join(outbase, "Annot_IS_SUCCESS"))):
        pathlib.Path(outbase).mkdir(exist_ok=True, parents=True)
        set_kmer_back =  defaultdict(set, {})
        dict_pept_backgrd = {}
        time_per_gene = [np.nan]
        mem_per_gene = [np.nan]
        all_gene_idxs = []

        for i, gene in enumerate(genes):
            ### measure time
            start_time = timeit.default_timer()

            # Genes not contained in the annotation...
            if (gene.name not in genetable.gene_to_cds_begin or \
                    gene.name not in genetable.gene_to_ts):
#                    logging.warning('>Gene {} is not in the genetable and not processed, please check the annotation file.'.format(gene.name))
                continue

            idx = get_idx(countinfo_b, output_sample, gene_idxs[i])
            all_gene_idxs.append(gene_idxs[i])

            chrm = gene.chr.strip()
            sub_mutation = get_sub_mutation_tuple(mutation, mutation_sample, chrm)
            ref_mut_seq = collect_background_transcripts(gene=gene, ref_seq_file=arg.ref_path, chrm=chrm, mutation=sub_mutation)

            # Gene counts information
            if arg.cross_graph_expr:
                countinfo_b = None
            if countinfo_b:
                gidx = countinfo_b.gene_idx_dict[gene.name]
                count_segments = np.arange(countinfo_b.gene_id_to_segrange[gidx][0], countinfo_b.gene_id_to_segrange[gidx][1])
                with h5py.File(countinfo_b.h5fname, 'r') as h5f:
                    seg_counts = h5f['segments'][count_segments, idx.sample]
            else:
                seg_counts = None

            get_and_write_background_peptide_and_kmer(peptide_dict = dict_pept_backgrd,
                                                      kmer_dict = set_kmer_back,
                                                      gene=gene,
                                                      ref_mut_seq=ref_mut_seq,
                                                       gene_table=genetable,
                                                       countinfo=countinfo_b,
                                                       seg_counts=seg_counts,
                                                       Idx=idx,
                                                       kmer=arg.kmer,
                                                       all_read_frames=arg.all_read_frames)

            time_per_gene.append(timeit.default_timer() - start_time)
            mem_per_gene.append(print_memory_diags(disable_print=True))

        save_backgrd_pep_dict(dict_pept_backgrd, filepointer, compression, outbase, verbose)
        dict_pept_backgrd.clear()
        for kmer_length in set_kmer_back:
            save_backgrd_kmer_set(set_kmer_back[kmer_length], filepointer, kmer_length, compression, outbase, verbose)
        set_kmer_back.clear()

        pathlib.Path(os.path.join(outbase, "Annot_IS_SUCCESS")).touch()
        logging_string = "....{}: annotation graph from batch {}/{} processed, max time cost: {}, memory cost:{} GB".format(
                                                                                output_sample,
                                                                                batch_name,
                                                                                len(gene_id_list),
                                                                                np.round(np.max(time_per_gene), 2),
                                                                                np.round(np.max(mem_per_gene), 2))
        logging.debug(logging_string)
        if (batch_name != 'all') and (batch_name % 10000 == 0):
            logging.info(logging_string)

    else:
        logging_string = "> {} : Batch {} exists, skip processing ".format(output_sample, batch_name)
        logging.debug(logging_string)
        if (batch_name != 'all') and (batch_name % 10000 == 0):
            logging.info(logging_string)

    return 'multiprocessing is success'



def process_gene_batch_foreground(output_sample, mutation_sample, output_samples_ids, genes, genes_info, gene_idxs, total_genes, genes_interest, disable_process_libsize, all_read_frames, complexity_cap, mutation, junction_dict, countinfo, genetable, arg, outbase, filepointer, compression, verbose):
    if arg.parallel > 1:
        batch_name = int(outbase.split('/')[-1].split('_')[-1])
    else:
        batch_name = 'all'
    if (arg.parallel==1) or (not os.path.exists(os.path.join(outbase, "output_sample_IS_SUCCESS"))):
        pathlib.Path(outbase).mkdir(exist_ok=True, parents=True)
        dict_kmer_foregr = defaultdict(dict, {})
        dict_pept_forgrd = {}
        time_per_gene = [np.nan]
        mem_per_gene = [np.nan]
        all_gene_idxs = []
        gene_expr = []
        
        for i, gene in enumerate(genes):
            ### measure time
            start_time = timeit.default_timer()

            # Genes not contained in the annotation in annotated CDS mode
            if (gene.name not in genetable.gene_to_cds_begin or \
                    gene.name not in genetable.gene_to_ts):
#                    logging.warning('>Gene {} is not in the genetable and not processed, please check the annotation file.'.format(gene.name))
                continue


            idx = get_idx(countinfo, output_sample, gene_idxs[i])
            # Gene counts information
            # Gene of interest always compute expression, others compute expession if required for library
            if not disable_process_libsize or (gene.name in genes_interest):
                if countinfo:
                    gidx = countinfo.gene_idx_dict[gene.name]

                    with h5py.File(countinfo.h5fname, 'r') as h5f:
                            # Get edges counts
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
                            # Get segment counts
                            seg_gene_idxs = np.arange(countinfo.gene_id_to_segrange[gidx][0],
                                                      countinfo.gene_id_to_segrange[gidx][1])
                            if arg.cross_graph_expr:
                                seg_counts = h5f['segments'][seg_gene_idxs, :]
                                if output_samples_ids is not None:
                                    seg_counts = seg_counts[:, output_samples_ids] # limitation fancy hdf5 indexing
                                else:
                                    output_samples_ids = np.arange(seg_counts.shape[1])
                            else:
                                seg_counts = h5f['segments'][seg_gene_idxs, idx.sample]
                else:
                    edge_idxs = None
                    edge_counts = None
                    seg_counts = None

            # library size calculated only for genes with CDS
            if (gene.name in genetable.gene_to_cds_begin) and (gene.name in genetable.gene_to_ts) and countinfo:
                gene_expr.append([gene.name] + get_total_gene_expr(gene, countinfo, idx, seg_counts, arg.cross_graph_expr))
            # Gene of interest
            if (gene.name not in genes_interest):
                continue
            # Genes with highly complex splicegraphs
            if (len(gene.splicegraph.vertices[1]) > complexity_cap):
                logging.warning('> Gene {} has a edge complexity > {}, not processed'.format(gene.name, complexity_cap))
                continue

            chrm = gene.chr.strip()
            sub_mutation = get_sub_mutation_tuple(mutation, mutation_sample, chrm)
            if (arg.mutation_sample is not None) and (arg.cross_graph_expr):
                mut_count_id = [idx for idx, sample in enumerate(arg.output_samples) if arg.mutation_sample.replace('-', '').replace('_', '').replace('.', '').replace('/', '') == sample][0]
            else:
                mut_count_id = None
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
                                            mut_count_id=mut_count_id,
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
                                            all_read_frames=arg.all_read_frames,
                                            filepointer=filepointer,
                                            graph_output_samples_ids = output_samples_ids,
                                            graph_samples=arg.output_samples,
                                            verbose_save=verbose
            )


            time_per_gene.append(timeit.default_timer() - start_time)
            mem_per_gene.append(print_memory_diags(disable_print=True))
            all_gene_idxs.append(gene_idxs[i])

        save_gene_expr_distr(gene_expr, arg.output_samples, output_sample,  filepointer, outbase, compression, verbose)
        save_forgrd_pep_dict(dict_pept_forgrd, filepointer, compression, outbase, arg.output_fasta, verbose)
        dict_pept_forgrd.clear()
        if not arg.cross_graph_expr:
            for kmer_length in dict_kmer_foregr:
                save_forgrd_kmer_dict(dict_kmer_foregr[kmer_length], filepointer, kmer_length, compression, outbase, verbose)
            dict_kmer_foregr.clear()
        if arg.cross_graph_expr and filepointer.kmer_segm_expr_fp['pqwriter'] is not None:
            filepointer.kmer_segm_expr_fp['pqwriter'].close()
            filepointer.kmer_edge_expr_fp['pqwriter'].close()

        pathlib.Path(os.path.join(outbase, "output_sample_IS_SUCCESS")).touch()
        logging_string = "....{}: output_sample graph from batch {}/{} processed, max time cost: {}, memory cost:{} GB".format(output_sample,
                                                                                          batch_name,
                                                                                          len(gene_id_list),
                                                                                          np.round(np.max(time_per_gene), 2),
                                                                                          np.round(np.max(mem_per_gene), 2))
        logging.debug(logging_string)
        if (batch_name != 'all') and (batch_name % 10000 == 0):
            logging.info(logging_string)

    else:
        logging_string = "> {} : Batch {} exists, skip processing ".format(output_sample, batch_name)
        logging.debug(logging_string)
        if (batch_name != 'all') and (batch_name % 10000 == 0):
            logging.info(logging_string)

    return 'multiprocessing is success'




def mode_build(arg):

    # read and process the annotation file
    logging.info(">>>>>>>>> Build: Start Preprocessing")
    logging.info('Building lookup structure ...')
    start_time = timeit.default_timer()
    genetable,chromosome_set = preprocess_ann(arg.ann_path)
    end_time = timeit.default_timer()
    logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()

    # load graph metadata
    start_time = timeit.default_timer()
    if arg.count_path is not None:
        logging.info('Loading count data ...')
        countinfo, matching_count_samples, matching_count_ids  = parse_gene_metadata_info(arg.count_path, arg.output_samples, arg.cross_graph_expr)

        end_time = timeit.default_timer()
        logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
        print_memory_diags()
        #size_factor = get_size_factor(strains, arg.libsize_path)
        #size_factor = None
    else:
        countinfo = None
        size_factor = None

    # read the variant file
    mutation = get_mutation_mode_from_parser(arg, matching_count_samples)

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
    #graph_data = graph_data[0:20]
    if arg.start_id != 0 and arg.start_id < len(graph_data): 
        logging.info("development feature: starting at gene number {}".format(arg.start_id))
        graph_data = graph_data[arg.start_id:]
        
    check_chr_consistence(chromosome_set,mutation,graph_data)


    # read the intron of interest file gtex_junctions.hdf5
    junction_dict = parse_junction_meta_info(arg.gtex_junction_path)

    # process the genes according to the annotation file
    # add CDS starts and reading frames to the respective nodes
    logging.info('Add reading frame to splicegraph ...')
    start_time = timeit.default_timer()
    graph_info, graph_data = genes_preprocess_all(graph_data, genetable.gene_to_cds_begin, arg.parallel, arg.all_read_frames)
    end_time = timeit.default_timer()
    logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()
    logging.info(">>>>>>>>> Finish Preprocessing")
    expr_distr_dict = {}

    # Parse user choice for genes
    graph_data, genes_interest, num, complexity_cap, disable_process_libsize = parse_gene_choices(arg.genes_interest, arg.process_chr, arg.process_num, arg.complexity_cap, arg.disable_process_libsize, graph_data)

    # process graph for each input output_sample
    output_libszie_fp = os.path.join(arg.output_dir,'expression_counts.libsize.tsv')
    logging.info(">>>>>>>>> Start traversing splicegraph")

    global gene_id_list
    global output_sample
    global filepointer
    
    # handle output_sample relatively to output mode 
    if arg.cross_graph_expr and countinfo:
        process_output_samples = ['cohort']
        # If samples requested, look for sample ids in count file
        if arg.output_samples:
            arg.output_samples = np.array(arg.output_samples)[np.argsort(matching_count_ids)]
            arg.output_samples = [output_sample.replace('-', '').replace('_', '').replace('.', '').replace('/', '')
                                  for output_sample in  arg.output_samples]
            output_samples_ids = matching_count_ids[np.argsort(matching_count_ids)]
        # If no samples requested, take all samples in countfile
        else:
            arg.output_samples = matching_count_samples
            output_samples_ids = None
    else:
        process_output_samples = arg.output_samples
        output_samples_ids = None
        
        
    for output_sample in process_output_samples:
        logging.info(">>>> Processing output_sample {}, there are {} graphs in total".format(output_sample,num))

        # prepare the output files
        if output_sample != arg.mutation_sample:
            output_path = os.path.join(arg.output_dir, '{}_mut{}'.format(output_sample, arg.mutation_sample))
        elif (output_sample == arg.mutation_sample) or (not arg.mutation_sample):
            output_path = os.path.join(arg.output_dir, output_sample)
        logging.info("Saving results to {}".format(output_path))

        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        gzip_tag = ''
        if arg.compressed:
            pq_compression = 'SNAPPY'
        else:
            pq_compression = None
        filepointer = initialize_fp(output_path, mutation.mode, gzip_tag,
                  arg.kmer, arg.output_fasta, arg.cross_graph_expr)

        # go over each gene in splicegraph
        gene_id_list = list(range(0,num))

        if arg.parallel > 1:

            logging.info('Parallel: {} Threads'.format(arg.parallel))
            batch_size = min(num, arg.batch_size)
            verbose_save = False
            gene_batches = [(i, gene_id_list[i:min(i + batch_size, len(gene_id_list))]) for i in
                            range(0, len(gene_id_list), batch_size)]

            if not arg.skip_annotation:
                # Build the background
                logging.info(">>>>>>>>> Start Background processing")
                pool_f = MyPool(processes=arg.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
                args = [(output_sample, arg.mutation_sample,  graph_data[gene_idx], gene_idx, arg.all_read_frames, mutation, countinfo, genetable, arg,
                      os.path.join(output_path, 'tmp_out_{}_batch_{}'.format(arg.mutation_mode, i + arg.start_id)), filepointer, None, verbose_save) for i, gene_idx in gene_batches ]

                result = pool_f.submit(mapper_funct_back, args)
                pool_f.terminate()
                exits_if_exception = [res for res in result]

            # Build the foreground
            logging.info(">>>>>>>>> Start Foreground processing")
            pool_f = MyPool(processes=arg.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
            args = [(output_sample, arg.mutation_sample, output_samples_ids, graph_data[gene_idx], graph_info[gene_idx], gene_idx, len(
                gene_id_list), genes_interest, disable_process_libsize, arg.all_read_frames, complexity_cap, mutation, junction_dict, countinfo, genetable, arg,
                  os.path.join(output_path, 'tmp_out_{}_batch_{}'.format(arg.mutation_mode, i + arg.start_id)), filepointer, None, verbose_save) for i, gene_idx in gene_batches ]

            result = pool_f.submit(mapper_funct, args)
            pool_f.terminate()
            exits_if_exception = [res for res in result]


            # Collects and pools the files of each batch
            logging.info("Start collecting results")
            if countinfo:
                collect_results(filepointer.gene_expr_fp, output_path, pq_compression, arg.mutation_mode)
            if arg.output_fasta:
                collect_results(filepointer.junction_peptide_fp, output_path, pq_compression, arg.mutation_mode)
            collect_results(filepointer.background_peptide_fp, output_path, pq_compression, arg.mutation_mode)
            collect_results(filepointer.junction_meta_fp, output_path, pq_compression, arg.mutation_mode)
            collect_results(filepointer.junction_kmer_fp, output_path, pq_compression, arg.mutation_mode, arg.kmer)
            collect_results(filepointer.background_kmer_fp, output_path, pq_compression, arg.mutation_mode, arg.kmer)
            collect_results(filepointer.kmer_segm_expr_fp, output_path, pq_compression, arg.mutation_mode)
            collect_results(filepointer.kmer_edge_expr_fp, output_path, pq_compression, arg.mutation_mode)
            if not arg.skip_tmpfiles_rm:
                logging.info("Cleaning temporary files")
                remove_folder_list(os.path.join(output_path, 'tmp_out_{}_batch'.format(arg.mutation_mode)))

        else:
            logging.info('Not Parallel')
            # Build the background
            logging.info(">>>>>>>>> Start Background processing")
            if not arg.skip_annotation:
                process_gene_batch_background(output_sample, arg.mutation_sample, graph_data, gene_id_list, arg.all_read_frames, mutation, countinfo, genetable, arg, output_path, filepointer, pq_compression, verbose=True)
            # Build the foreground and remove the background if needed
            logging.info(">>>>>>>>> Start Foreground processing")
            process_gene_batch_foreground( output_sample, arg.mutation_sample, output_samples_ids, graph_data, graph_info, gene_id_list, len(gene_id_list), genes_interest, disable_process_libsize, arg.all_read_frames, complexity_cap, mutation, junction_dict,
                             countinfo, genetable, arg, output_path, filepointer, pq_compression, verbose=True)

        if (not disable_process_libsize) and countinfo:
            create_libsize(filepointer.gene_expr_fp, output_libszie_fp, output_sample)




