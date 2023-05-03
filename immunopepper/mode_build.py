# Python libraries
""""""
# Core operation of ImmunoPepper. Traverse splicegraph and get kmer/peptide output
from collections import defaultdict
import h5py
import logging
import multiprocessing as mp
from multiprocessing.pool import Pool
from multiprocessing.pool import ThreadPool
import numpy as np
import os
import pathlib
import pickle
import signal as sig
import sys
import timeit

# immuno module
from immunopepper.io_ import collect_results
from immunopepper.io_ import get_save_path
from immunopepper.io_ import initialize_fp
from immunopepper.io_ import remove_folder_list
from immunopepper.io_ import save_bg_kmer_set
from immunopepper.io_ import save_bg_peptide_set
from immunopepper.io_ import save_gene_expr_distr
from immunopepper.io_ import save_fg_kmer_dict
from immunopepper.io_ import save_fg_peptide_set

from immunopepper.mutations import get_sub_mutations
from immunopepper.mutations import load_mutations
from immunopepper.preprocess import genes_preprocess_all
from immunopepper.preprocess import parse_junction_meta_info
from immunopepper.preprocess import parse_gene_choices
from immunopepper.preprocess import parse_gene_metadata_info
from immunopepper.preprocess import parse_output_samples_choices
from immunopepper.preprocess import parse_uniprot
from immunopepper.preprocess import preprocess_ann
from immunopepper.traversal import collect_background_transcripts
from immunopepper.traversal import collect_vertex_pairs
from immunopepper.traversal import get_and_write_background_peptide_and_kmer
from immunopepper.traversal import get_and_write_peptide_and_kmer
from immunopepper.utils import check_chr_consistence
from immunopepper.utils import create_libsize
from immunopepper.utils import get_idx
from immunopepper.utils import get_total_gene_expr
from immunopepper.utils import print_memory_diags


### intermediate fix to load pickle files stored under previous version
from spladder.classes import gene as cgene
from spladder.classes import splicegraph as csplicegraph
from spladder.classes import segmentgraph as csegmentgraph
sys.modules['modules'] = cgene
sys.modules['modules.classes.gene'] = cgene
sys.modules['modules.classes.splicegraph'] = csplicegraph
sys.modules['modules.classes.segmentgraph'] = csegmentgraph
### end fix

def pool_initializer_glob(countinfo_glob, genetable_glob, kmer_database_glob): #Moved from utils because of global variables
    global countinfo
    global genetable
    global kmer_database
    countinfo = countinfo_glob
    genetable = genetable_glob
    kmer_database = kmer_database_glob
    #return sig.signal(sig.SIGINT, sig.SIG_IGN)

def mapper_funct(tuple_arg):
    process_gene_batch_foreground(*tuple_arg)


def mapper_funct_back(tuple_arg):
    process_gene_batch_background(*tuple_arg)


def process_gene_batch_background(output_sample, mutation_sample, genes, gene_idxs, n_genes, mutation,
                                  genetable, arg, outbase, filepointer, verbose=False):
    if arg.parallel > 1:
        batch_name = int(outbase.split('/')[-1].split('_')[-1])
    else:
        batch_name = 'all'

    if (arg.parallel==1) or (not os.path.exists(os.path.join(outbase, "Annot_IS_SUCCESS"))):
        pathlib.Path(outbase).mkdir(exist_ok=True, parents=True)
        set_pept_backgrd = set()
        set_kmer_backgrd = set()
        time_per_gene = []
        mem_per_gene = []
        all_gene_idxs = []

        for i, gene in enumerate(genes):
            ### measure time
            start_time = timeit.default_timer()

            # Genes not contained in the annotation...
            if (gene.name not in genetable.gene_to_cds_begin or \
                    gene.name not in genetable.gene_to_ts):
                continue

            all_gene_idxs.append(gene_idxs[i])

            chrm = gene.chr.strip()
            sub_mutation = get_sub_mutations(mutation, mutation_sample, chrm)
            ref_mut_seq = collect_background_transcripts(gene=gene, ref_seq_file=arg.ref_path,
                                                         chrm=chrm, mutation=sub_mutation)


            get_and_write_background_peptide_and_kmer(peptide_set=set_pept_backgrd,
                                                      kmer_set=set_kmer_backgrd,
                                                      gene=gene,
                                                      ref_mut_seq=ref_mut_seq,
                                                      gene_table=genetable,
                                                      countinfo=None,
                                                      kmer_length=arg.kmer,
                                                      all_read_frames=arg.all_read_frames)

            time_per_gene.append(timeit.default_timer() - start_time)
            mem_per_gene.append(print_memory_diags(disable_print=True))

        save_bg_peptide_set(set_pept_backgrd, filepointer, outbase, verbose)
        set_pept_backgrd.clear()
        save_bg_kmer_set(set_kmer_backgrd, filepointer, outbase, verbose)
        set_kmer_backgrd.clear()

        pathlib.Path(os.path.join(outbase, "Annot_IS_SUCCESS")).touch()

        if time_per_gene:
            logging_string = (f'....{output_sample}: annotation graph from batch {batch_name}/{n_genes} '
                              f'processed, max time cost: {np.round(np.nanmax(time_per_gene), 2)}, '
                              f'memory cost: {np.round(np.nanmax(mem_per_gene), 2)} GB')
            logging.debug(logging_string)
        else:
            logging_string = (f'....{output_sample}: output_sample graph from batch {batch_name}/{n_genes}, no processing')

        if (batch_name != 'all') and (batch_name % 10000 == 0):
            logging.info(logging_string)

    else:  # Batch has already been saved to disk
        logging_string = f'> {output_sample} : Batch {batch_name} exists, skip processing'
        logging.debug(logging_string)
        if (batch_name != 'all') and (batch_name % 10000 == 0):
            logging.info(logging_string)

    return 'multiprocessing is success'


def process_gene_batch_foreground(output_sample, mutation_sample, output_samples_ids, genes,
                                  genes_info, gene_idxs, n_genes, genes_interest, disable_process_libsize,
                                  all_read_frames, complexity_cap, mutation, junction_dict,
                                  arg, outbase, filepointer, verbose):
    global countinfo
    global genetable
    global kmer_database
    if arg.parallel > 1:
        batch_name = int(outbase.split('/')[-1].split('_')[-1])
    else:
        batch_name = 'all'

    if (arg.parallel==1) or (not os.path.exists(os.path.join(outbase, "output_sample_IS_SUCCESS"))):
        pathlib.Path(outbase).mkdir(exist_ok=True, parents=True)
        set_pept_forgrd = set()
        time_per_gene = []
        mem_per_gene = []
        all_gene_idxs = []
        gene_expr = []

        for i, gene in enumerate(genes):
            ### measure time
            start_time = timeit.default_timer()

            # Genes not contained in the annotation in annotated CDS mode
            if (gene.name not in genetable.gene_to_cds_begin or \
                    gene.name not in genetable.gene_to_ts):
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
                                edge_gene_idxs = np.arange(countinfo.gene_id_to_edgerange[gidx][0],
                                                           countinfo.gene_id_to_edgerange[gidx][1])
                                edge_idxs = h5f['edge_idx'][list(edge_gene_idxs)].astype('int')
                                edge_counts = h5f['edges'][edge_gene_idxs,:] # will compute expression on whole graph

                            # Get segment counts
                            seg_gene_idxs = np.arange(countinfo.gene_id_to_segrange[gidx][0],
                                                      countinfo.gene_id_to_segrange[gidx][1])
                            seg_counts = h5f['segments'][seg_gene_idxs, :]
                            if output_samples_ids is not None:
                                seg_counts = seg_counts[:, output_samples_ids] # limitation fancy hdf5 indexing
                            else:
                                output_samples_ids = np.arange(seg_counts.shape[1])
                else:
                    edge_idxs = None
                    edge_counts = None
                    seg_counts = None

            # library size calculated only for genes with CDS
            if (gene.name in genetable.gene_to_cds_begin) and (gene.name in genetable.gene_to_ts) and countinfo:
                gene_expr.append([gene.name] + get_total_gene_expr(gene, countinfo, seg_counts))

            # Process only gene quantifications and libsizes
            if arg.libsize_extract:
                time_per_gene.append(timeit.default_timer() - start_time)
                mem_per_gene.append(print_memory_diags(disable_print=True))
                continue

            # Process only gene of interest
            if (gene.name not in genes_interest):
                continue

            # Do not process genes with highly complex splicegraphs
            if (len(gene.splicegraph.vertices[1]) > complexity_cap):
                logging.warning(f'> Gene {gene.name} has a edge complexity > {complexity_cap}, not processed')
                continue

            chrm = gene.chr.strip()
            sub_mutation = get_sub_mutations(mutation, mutation_sample, chrm)
            if (arg.mutation_sample is not None):
                mut_count_id = [idx for idx, sample in enumerate(arg.output_samples)
                                if arg.mutation_sample.replace('-', '').replace('_', '').replace('.', '').replace('/', '') == sample][0]
            else:
                mut_count_id = None
            junction_list = None
            if not junction_dict is None and chrm in junction_dict:
                junction_list = junction_dict[chrm]
 
            pathlib.Path(get_save_path(filepointer.kmer_segm_expr_fp, outbase)).mkdir(exist_ok=True, parents=True)
            pathlib.Path(get_save_path(filepointer.kmer_edge_expr_fp, outbase)).mkdir(exist_ok=True, parents=True)
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
                                                 kmer_length=arg.kmer,
                                                 filter_redundant=arg.filter_redundant)

            get_and_write_peptide_and_kmer(peptide_set=set_pept_forgrd,
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
                                            kmer_database=kmer_database,
                                            kmer=arg.kmer,
                                            force_ref_peptides=arg.force_ref_peptides,
                                            out_dir=outbase,
                                            edge_idxs=edge_idxs,
                                            edge_counts=edge_counts,
                                            seg_counts=seg_counts,
                                            all_read_frames=arg.all_read_frames,
                                            filepointer=filepointer,
                                            graph_output_samples_ids = output_samples_ids,
                                            graph_samples=arg.output_samples,
                                            verbose_save=verbose
            )


            time_per_gene.append(timeit.default_timer() - start_time)
            mem_per_gene.append(print_memory_diags(disable_print=True))
            all_gene_idxs.append(gene_idxs[i])

        save_gene_expr_distr(gene_expr, arg.output_samples, output_sample,  filepointer, outbase, verbose)
        save_fg_peptide_set(set_pept_forgrd, filepointer, outbase, arg.output_fasta, verbose)
        set_pept_forgrd.clear()

        pathlib.Path(os.path.join(outbase, "output_sample_IS_SUCCESS")).touch()

        if time_per_gene:
            logging_string = (f'....{output_sample}: output_sample graph from batch {batch_name}/{n_genes} processed, '
                          f'max time cost: {np.round(np.nanmax(time_per_gene), 2)}, '
                          f'memory cost: {np.round(np.nanmax(mem_per_gene), 2)} GB')
            logging.debug(logging_string)
        else:
            logging_string = (f'....{output_sample}: output_sample graph from batch {batch_name}/{n_genes}, no processing')
        if (batch_name != 'all') and (batch_name % 10000 == 0):
            logging.info(logging_string)

    else: # Batch has already been saved to disk
        logging_string = f'> {output_sample}: Batch {batch_name} exists, skip processing'
        logging.debug(logging_string)
        if (batch_name != 'all') and (batch_name % 10000 == 0):
            logging.info(logging_string)

    return 'multiprocessing is success'


def mode_build(arg):
    global output_sample
    global filepointer
    global countinfo #Will be used in non parallel mode
    global genetable #Will be used in non parallel mode
    global kmer_database #Will be used in non parallel mode
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
        countinfo, matching_count_samples, matching_count_ids  = parse_gene_metadata_info(arg.count_path,
                                                                                          arg.output_samples)


        end_time = timeit.default_timer()
        logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
        print_memory_diags()
        #size_factor = get_size_factor(strains, arg.libsize_path)
    else:
        countinfo = None
        matching_count_samples = None
        matching_count_ids = None

    # read the variant file

    mutation = load_mutations(arg.germline, arg.somatic, arg.mutation_sample, arg.heter_code, 
                              arg.pickle_samples if arg.use_mut_pickle else None,
                              arg.sample_name_map, arg.output_dir if arg.use_mut_pickle else None)

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
    #graph_data = graph_data[940:942]
    if arg.start_id != 0 and arg.start_id < len(graph_data):
        logging.info(f'development feature: starting at gene number {arg.start_id}')
        graph_data = graph_data[arg.start_id:]

    check_chr_consistence(chromosome_set, mutation, graph_data)


    # read the intron of interest file gtex_junctions.hdf5
    junction_dict = parse_junction_meta_info(arg.gtex_junction_path)

    # read and process uniprot file
    kmer_database = parse_uniprot(arg.kmer_database)

    # process the genes according to the annotation file
    # add CDS starts and reading frames to the respective nodes
    logging.info('Add reading frame to splicegraph ...')
    start_time = timeit.default_timer()
    graph_info, graph_data = genes_preprocess_all(graph_data, genetable.gene_to_cds_begin,
                                                  arg.parallel, arg.all_read_frames)
    end_time = timeit.default_timer()
    logging.info('\tTime spent: {:.3f} seconds'.format(end_time - start_time))
    print_memory_diags()
    logging.info(">>>>>>>>> Finish Preprocessing")
    expr_distr_dict = {}

    # parse user choice for genes
    graph_data, genes_interest, n_genes, \
    complexity_cap, disable_process_libsize = parse_gene_choices(arg.genes_interest, arg.process_chr, arg.process_num,
                                                                 arg.complexity_cap, arg.disable_process_libsize,
                                                                 graph_data)

    # process graph for each input output_sample
    output_libsize_fp = os.path.join(arg.output_dir, 'expression_counts.libsize.tsv')


    # parse output_sample relatively to output mode
    process_output_samples, output_samples_ids = parse_output_samples_choices(arg, countinfo, matching_count_ids,
                                                                              matching_count_samples)
    logging.info(">>>>>>>>> Start traversing splicegraph")
    for output_sample in process_output_samples:
        logging.info(f'>>>> Processing output_sample {output_sample}, there are {n_genes} graphs in total')

        # prepare the output files
        if output_sample != arg.mutation_sample:
            output_path = os.path.join(arg.output_dir, f'{output_sample}_mut{arg.mutation_sample}')
        elif (output_sample == arg.mutation_sample) or (not arg.mutation_sample):
            output_path = os.path.join(arg.output_dir, output_sample)
        logging.info(f'Saving results to {output_path}')

        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        gzip_tag = '' #TODO clean next commit
        if arg.compressed:
            pq_compression = 'SNAPPY'
        else:
            pq_compression = None
        filepointer = initialize_fp(output_path, mutation.mode, arg.output_fasta)

        # go over each gene in splicegraph
        genes_range = list(range(0, n_genes))

        if arg.parallel > 1:
            logging.info(f'Parallel: {arg.parallel} Threads')
            batch_size = min(n_genes, arg.batch_size)
            verbose_save = False
            gene_batches = [(i, genes_range[i:min(i + batch_size, n_genes)]) for i in
                            range(0, n_genes, batch_size)]

            if (not arg.skip_annotation) and not (arg.libsize_extract):
                # Build the background
                logging.info(">>>>>>>>> Start Background processing")
                with ThreadPool(processes=None, initializer=pool_initializer_glob, initargs=(countinfo, genetable, kmer_database)) as pool:
                    args = [(output_sample, arg.mutation_sample,  graph_data[gene_idx], gene_idx, n_genes, mutation,
                             genetable, arg,
                             os.path.join(output_path, f'tmp_out_{mutation.mode}_batch_{i + arg.start_id}'),
                             filepointer, verbose_save) for i, gene_idx in gene_batches ]
                    result = pool.imap(mapper_funct_back, args, chunksize=1)
                    exits_if_exception = [res for res in result]

            # Build the foreground
            logging.info(">>>>>>>>> Start Foreground processing")
            with Pool(processes=None, initializer=pool_initializer_glob, initargs=(countinfo, genetable, kmer_database)) as pool:
                args = [(output_sample, arg.mutation_sample, output_samples_ids, graph_data[gene_idx],
                         graph_info[gene_idx], gene_idx, n_genes, genes_interest, disable_process_libsize,
                         arg.all_read_frames, complexity_cap, mutation, junction_dict, arg,
                         os.path.join(output_path, f'tmp_out_{mutation.mode}_batch_{i + arg.start_id}'),
                         filepointer, verbose_save) for i, gene_idx in gene_batches ]
                result = pool.imap(mapper_funct, args, chunksize=1)
                exits_if_exception = [res for res in result]


            logging.info("Finished traversal")


        else:
            logging.info('Not Parallel')
            # Build the background
            logging.info(">>>>>>>>> Start Background processing")
            if (not arg.skip_annotation) and not (arg.libsize_extract):
                process_gene_batch_background(output_sample, arg.mutation_sample, graph_data, genes_range, n_genes,
                                              mutation, genetable, arg, output_path, filepointer,
                                              verbose=True)
            # Build the foreground and remove the background if needed
            logging.info(">>>>>>>>> Start Foreground processing")
            process_gene_batch_foreground( output_sample, arg.mutation_sample, output_samples_ids, graph_data,
                                           graph_info, genes_range, n_genes, genes_interest, disable_process_libsize,
                                           arg.all_read_frames, complexity_cap, mutation, junction_dict,
                                           arg, output_path, filepointer,
                                           verbose=True)

        if (not disable_process_libsize) and countinfo:
            create_libsize(filepointer.gene_expr_fp, output_libsize_fp, output_path, mutation.mode, arg.parallel)




