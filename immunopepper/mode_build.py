# Python libraries
""""""
# Core operation of ImmunoPepper. Traverse splicegraph and get kmer/peptide output
from collections import defaultdict
import h5py
import logging
import multiprocessing as mp
import numpy as np
import os
import pathlib
import pickle
import signal as sig
import sys
import timeit

# immuno module
from .constant import NOT_EXIST
from .filter import add_set_kmer_back
from .filter import add_dict_kmer_forgrd
from .filter import add_dict_peptide
from .io_ import collect_results
from .io_ import gz_and_normal_open
from .io_ import initialize_fp
from .io_ import remove_folder_list
from .io_ import save_backgrd_kmer_set
from .io_ import save_backgrd_pep_dict
from .io_ import save_forgrd_kmer_dict
from .io_ import save_forgrd_pep_dict
from .io_ import write_gene_expr
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


def process_gene_batch_background(sample, genes, gene_idxs,  mutation , countinfo, genetable, arg, outbase, remove_annot, uniq_foreground, compression=None, verbose=False):
    results = []
    for i, gene in enumerate(genes):
        ### measure time
        start_time = timeit.default_timer()

        ### set result dict
        R = dict()
        R['gene_name'] = gene.name

        # Genes not contained in the annotation...
        if gene.name not in genetable.gene_to_cds_begin or \
                gene.name not in genetable.gene_to_ts:
            #logger.warning('>Gene name {} is not in the genetable and not processed, please check the annotation file.'.format(gene.name))
            R['processed'] = False
            results.append(R)
            continue
        R['processed'] = True

        idx = get_idx(countinfo, sample, gene_idxs[i])
        R['gene_idx'] = gene_idxs[i]

        chrm = gene.chr.strip()
        sub_mutation = get_sub_mutation_tuple(mutation, sample, chrm)
        ref_mut_seq = collect_background_transcripts(gene=gene, ref_seq_file=arg.ref_path, chrm=chrm, mutation=sub_mutation)

        # Gene counts information
        if countinfo:
            gidx = countinfo.gene_idx_dict[gene.name]
            count_segments = np.arange(countinfo.gene_id_to_segrange[gidx][0], countinfo.gene_id_to_segrange[gidx][1])
            with h5py.File(countinfo.h5fname, 'r') as h5f:
                seg_counts = h5f['segments'][count_segments, idx.sample]

        R['background_peptide_list'], \
        R['background_kmer_lists'] = get_and_write_background_peptide_and_kmer(gene=gene,
                                                                               ref_mut_seq=ref_mut_seq,
                                                                               gene_table=genetable,
                                                                               countinfo=countinfo,
                                                                               seg_counts=seg_counts,
                                                                               Idx=idx,
                                                                               kmer=arg.kmer)

        R['time'] = timeit.default_timer() - start_time
        R['memory'] = print_memory_diags(disable_print=True)
        R['outbase'] = outbase
        results.append(R)
    process_result(results, ['background_peptide_list', 'background_kmer_lists'], outbase, remove_annot, uniq_foreground, compression, verbose)
    return dict_pept_backgrd, set_kmer_back  #Does not update the dictionaries and list globally because of the subprocess


def process_gene_batch_foreground(sample, genes, genes_info, gene_idxs, total_genes, all_read_frames, mutation, junction_dict, countinfo, genetable, arg, outbase, dict_pept_backgrd, remove_annot, uniq_foreground, compression, verbose):
    results = []
    for i, gene in enumerate(genes):
        ### measure time
        start_time = timeit.default_timer()

        ### set result dict
        R = dict()
        R['gene_name'] = gene.name

        # Genes not contained in the annotation...
        if gene.name not in genetable.gene_to_cds_begin or \
           gene.name not in genetable.gene_to_ts:

            #logger.warning('>Gene name {} is not in the genetable and not processed, please check the annotation file.'.format(gene.name))
            R['processed'] = False
            results.append(R)
            continue
        R['processed'] = True


        idx = get_idx(countinfo, sample, gene_idxs[i])

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
                    edge_counts = h5f['edges'][edge_gene_idxs, idx.sample]
                seg_gene_idxs = np.arange(countinfo.gene_id_to_segrange[gidx][0],
                                          countinfo.gene_id_to_segrange[gidx][1])
                seg_counts = h5f['segments'][seg_gene_idxs, idx.sample]

        R['total_expr'] = get_total_gene_expr(gene, countinfo, idx, seg_counts)
        R['gene_idx'] = gene_idxs[i]

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


        R['output_metadata_list'], \
        R['output_kmer_lists'] = get_and_write_peptide_and_kmer(gene=gene,
                                                                vertex_pairs=vertex_pairs,
                                                                background_pep_dict=dict_pept_backgrd,
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
                                                                seg_counts=seg_counts
                                                                )
        # Gene expression and sample library size
        gene_name_expr_distr.append((R['gene_name'], R['total_expr']))
        expr_distr.append(R['total_expr'])

        R['time'] = timeit.default_timer() - start_time
        R['memory'] = print_memory_diags(disable_print=True)
        R['outbase'] = outbase
        results.append(R)

    process_result(results, ['output_metadata_list', 'output_kmer_lists'], outbase, remove_annot, uniq_foreground, compression, verbose)
    return gene_name_expr_distr, expr_distr,  dict_pept_forgrd, dict_kmer_foregr #Does not update the dictionaries and list globally because of the subprocess


def process_result(gene_results, output_name, outbase, remove_annot=False, uniq_foreground=False, compression=None, verbose=True):
    '''
    Parameters
    -----------
    output_name:
    output_name: a list which contains the output types. E.g. 'background_peptide_list', 'background_kmer_lists', 'output_metadata_list', 'output_kmer_lists'
    remove_annot: whether to remove the background

    '''
    global sample
    global dict_kmer_foregr
    global set_kmer_back
    global dict_pept_forgrd
    global dict_pept_backgrd
    global filepointer
    pool_genes = defaultdict(list, {})
    time_per_gene = []
    mem_per_gene = []
    gene_idxs = []
    for gene_result in gene_results:
        if gene_result['processed']:
            for result_type in output_name:
                pool_genes[result_type].extend(gene_result[result_type])
                time_per_gene.append(gene_result['time'])
                mem_per_gene.append(gene_result['memory'])
                gene_idxs.append(gene_result['gene_idx'])
    s1 = timeit.default_timer()
    write_gene_result(pool_genes, dict_pept_forgrd, dict_pept_backgrd, dict_kmer_foregr, set_kmer_back,
                      filepointer, remove_annot, uniq_foreground, compression, outbase, verbose)

    if gene_idxs:
        logging.info("> {}: {}/{} processed, mean time cost: {}, memory cost:{} GB for per gene; writing results took {} seconds ".format(sample, gene_idxs[-1]  + 1,
                                                                                  len(gene_id_list),
                                                                                  np.mean(time_per_gene),
                                                                                  np.mean(mem_per_gene), 
                                                                                  timeit.default_timer() - s1))
    del gene_results



def write_gene_result(gene_result, dict_pept_forgrd, dict_pept_backgrd, dict_kmer_foregr, set_kmer_back, filepointer,
                      remove_annot=True, uniq_foreground=False, compression=None, outbase=None,  verbose=True):

    if ('background_peptide_list' in gene_result) and (len(gene_result['background_peptide_list'])):
        records = gene_result['background_peptide_list']
        add_dict_peptide(dict_pept_backgrd, records)
        save_backgrd_pep_dict(dict_pept_backgrd, filepointer, compression, outbase, verbose)

    if ('background_kmer_lists' in gene_result) and (len(gene_result['background_kmer_lists'])):
        records = [item for sublist in gene_result['background_kmer_lists'] for item in sublist]
        add_set_kmer_back(set_kmer_back, records)
        if not remove_annot:
            save_backgrd_kmer_set(set_kmer_back, filepointer, compression, outbase, verbose)
            set_kmer_back.clear()


    if ('output_metadata_list' in gene_result) and (len(gene_result['output_metadata_list'])):
        records = gene_result['output_metadata_list']
        add_dict_peptide(dict_pept_forgrd, records)
        if not uniq_foreground:
            save_forgrd_pep_dict(dict_pept_forgrd, filepointer, compression, outbase, verbose)
            dict_pept_forgrd.clear()

    if ('output_kmer_lists' in gene_result) and (len(gene_result['output_kmer_lists'])):
        records = [item for sublist in gene_result['output_kmer_lists'] for item in sublist]
        add_dict_kmer_forgrd(dict_kmer_foregr, records,  set_kmer_back, remove_annot)
        if not uniq_foreground:
            save_forgrd_kmer_dict(dict_kmer_foregr, filepointer, compression, outbase, verbose)
            dict_kmer_foregr.clear()




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
    all_read_frames = arg.all_read_frames
    remove_annot =  arg.remove_annot
    uniq_foreground = arg.uniq_foreground or remove_annot
    if uniq_foreground:
        logging.info('INFO: kmers and peptides will be kept in memory')

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
    global gene_name_expr_distr
    global expr_distr
    global gene_id_list
    global sample
    global dict_kmer_foregr
    global set_kmer_back
    global dict_pept_forgrd
    global dict_pept_backgrd
    global filepointer


    for sample in arg.samples:
        logging.info(">>>> Processing sample {}, there are {} graphs in total".format(sample,num))
        error_gene_num = 0 #TODO update
        time_list = []
        memory_list = []
        expr_distr = []
        gene_name_expr_distr = []
        dict_kmer_foregr = {}
        set_kmer_back = set()
        dict_pept_forgrd = {}
        dict_pept_backgrd = {}

        # prepare for the output file
        output_path = os.path.join(arg.output_dir, sample)
        if not os.path.isdir(output_path):
            os.makedirs(output_path)

        junction_meta_file_path = os.path.join(output_path, mutation.mode + '_metadata.tsv.gz.pq')

        if arg.compressed:
            gzip_tag = '.gz'
            pq_compression = 'GZIP'
        else:
            pq_compression = None
            gzip_tag = ''
        junction_peptide_file_path = os.path.join(output_path, mutation.mode + '_peptides.fa.pq'+gzip_tag)
        background_peptide_file_path = os.path.join(output_path, mutation.mode + '_back_peptides.fa.pq'+gzip_tag)
        junction_kmer_file_path = os.path.join(output_path, mutation.mode + '_junction_kmer.pq'+gzip_tag)
        background_kmer_file_path = os.path.join(output_path, mutation.mode + '_back_kmer.pq'+gzip_tag)
        gene_expr_file_path = os.path.join(output_path, 'gene_expression_detail.tsv'+gzip_tag)
        gene_expr_fp = gz_and_normal_open(gene_expr_file_path, 'w')


        # go over each gene in splicegraph
        gene_id_list = list(range(0,num))

        if arg.parallel > 1:
            filepointer = initialize_fp(junction_peptide_file_path, junction_meta_file_path,
                                        background_peptide_file_path,
                                        junction_kmer_file_path, background_kmer_file_path,
                                        compression=pq_compression)
            logging.info('Parallel: {} Threads'.format(arg.parallel))


            batch_size = min(num, arg.batch_size)
            verbose_save = False
            # Build the background
            logging.info(">>>>>>>>> Start Background processing")
            pool = mp.Pool(processes=arg.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
            for i in range(0, len(gene_id_list), batch_size):
                gene_idx = gene_id_list[i:min(i + batch_size, len(gene_id_list))]
                outbase = os.path.join(output_path, 'tmp_out_{}_{}'.format(arg.mutation_mode, i))
                pathlib.Path(outbase).mkdir(exist_ok= True, parents= True)
                res = pool.apply_async(process_gene_batch_background, args=(sample, graph_data[gene_idx], gene_idx, mutation, countinfo, genetable, arg, outbase, remove_annot, uniq_foreground, pq_compression, verbose_save))
            pool.close()
            pool.join()
            dict_pept_backgrd, set_kmer_back = res.get()

            # Build the foreground and remove the background if needed
            logging.info(">>>>>>>>> Start Foreground processing")
            pool = mp.Pool(processes=arg.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
            for i in range(0, len(gene_id_list), batch_size):
                gene_idx = gene_id_list[i:min(i + batch_size, len(gene_id_list))]
                outbase = os.path.join(output_path, 'tmp_out_{}_{}'.format(arg.mutation_mode, i))
                res = pool.apply_async(process_gene_batch_foreground, args=(sample, graph_data[gene_idx], graph_info[gene_idx], gene_idx, len(gene_id_list), all_read_frames, mutation, junction_dict, countinfo, genetable, arg, outbase, dict_pept_backgrd, remove_annot, uniq_foreground, pq_compression, verbose_save))
            pool.close()
            pool.join()
            gene_name_expr_distr, expr_distr, dict_pept_forgrd, dict_kmer_foregr = res.get()

            # Collects and pools the files of each batch
            logging.debug('start collecting results')
            collect_results(filepointer.background_peptide_fp, output_path, pq_compression)
            collect_results(filepointer.junction_peptide_fp,output_path, pq_compression)
            collect_results(filepointer.junction_meta_fp,output_path, pq_compression)
            collect_results(filepointer.junction_kmer_fp,output_path, pq_compression)
            collect_results(filepointer.background_kmer_fp,output_path, pq_compression)
            remove_folder_list(os.path.join(output_path, 'tmp_out_')) #TODO add back after development

        else:
            filepointer = initialize_fp(junction_peptide_file_path, junction_meta_file_path,
                                        background_peptide_file_path,
                                        junction_kmer_file_path, background_kmer_file_path,
                                        compression=pq_compression)
            logging.info('Not Parallel')
            # Build the background
            process_gene_batch_background(sample, graph_data, gene_id_list, mutation, countinfo, genetable, arg, output_path, remove_annot, uniq_foreground, pq_compression, verbose=True)
            # Build the foreground and remove the background if needed
            process_gene_batch_foreground( sample, graph_data, graph_info, gene_id_list, len(gene_id_list), all_read_frames, mutation, junction_dict,
                             countinfo, genetable, arg, output_path, dict_pept_backgrd, remove_annot, uniq_foreground, pq_compression, verbose=True)


        expr_distr_dict[sample] = expr_distr
        write_gene_expr(gene_expr_fp,gene_name_expr_distr)
        create_libsize(expr_distr_dict,output_libszie_fp)


        # Save the data structures kept in memory for filtering
        if uniq_foreground:
            save_forgrd_pep_dict(dict_pept_forgrd, filepointer, pq_compression, outbase=None, verbose=True)
            del dict_pept_forgrd

            save_forgrd_kmer_dict(dict_kmer_foregr, filepointer, pq_compression, outbase=None, verbose=True)
            del dict_kmer_foregr

        if remove_annot:
            save_backgrd_pep_dict(dict_pept_backgrd , filepointer, pq_compression, outbase=None, verbose=True)
            del dict_pept_backgrd

            save_backgrd_kmer_set(set_kmer_back, filepointer, pq_compression, outbase=None, verbose=True)
            del set_kmer_back


    if uniq_foreground:
        uniq_apply = "Foreground kmers/peptides are made unique across genes and the metadata is aggregated per kmer/peptide"
    else:
        uniq_apply = "Foreground kmers/peptides were not made unique across genes. Metadata is aggregated for kmer/peptide within each gene"
    if remove_annot:
        filt_apply = "Removed kmers present in annotation from the kmer foreground output"
    else:
        filt_apply = "Did not apply filtering for annotation on foreground kmers"
    logging.info(">>>>>>>>> Build: Finish traversing splicegraph in mutation mode {}')".format(mutation.mode))
    logging.info(">>>>>>>>> {}".format(filt_apply))
    logging.info(">>>>>>>>> {}".format(uniq_apply))

