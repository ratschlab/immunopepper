# Python libraries
""""""
# Core operation of ImmunoPepper. Traverse splicegraph and get kmer/peptide output
from collections import defaultdict
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
from .io_ import gz_and_normal_open
from .io_ import write_gene_expr
from .mutations import get_mutation_mode_from_parser
from .mutations import get_sub_mutation_tuple
from .preprocess import genes_preprocess_all
from .preprocess import parse_junction_meta_info
from .preprocess import parse_gene_metadata_info
from .preprocess import preprocess_ann
from .traversal import collect_vertex_pairs
from .traversal import get_and_write_background_peptide_and_kmer
from .traversal import get_and_write_peptide_and_kmer
from .inmemory import collect_results
from .inmemory import filter_onkey_dict
from .inmemory import initialize_parquet
from .inmemory import remove_folder_list
from .inmemory import save_backgrd_kmer_dict
from .inmemory import save_forgrd_kmer_dict
from .inmemory import save_forgrd_pep_dict
from .inmemory import write_gene_result
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


def process_gene_batch(sample, genes, genes_info, gene_idxs, total_genes, mutation, junction_dict, countinfo, genetable, arg, outbase):
    results = []
    for i, gene in enumerate(genes):
        ### measure time
        start_time = timeit.default_timer()

        ### set result dict
        R = dict()
        R['gene_name'] = gene.name

        # Genes not contained in the annotation...
        if gene.name not in genetable.gene_to_cds_begin or \
           gene.name not in genetable.gene_to_ts or \
           countinfo.gene_idx_dict[gene.name] not in countinfo.gene_id_to_edgerange or \
           countinfo.gene_idx_dict[gene.name] not in countinfo.gene_id_to_segrange:
            #logger.warning('>Gene name {} is not in the genetable and not processed, please check the annotation file.'.format(gene.name))
            R['processed'] = False
            results.append(R)
            continue
        R['processed'] = True

        idx = get_idx(countinfo, sample, gene_idxs[i])
        R['total_expr'] = get_total_gene_expr(gene, countinfo, idx)
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
                                             disable_concat=arg.disable_concat,
                                             kmer=arg.kmer,
                                             filter_redundant=arg.filter_redundant)

        R['background_peptide_list'], \
        R['background_kmer_lists'] = get_and_write_background_peptide_and_kmer(gene=gene, 
                                                                               ref_mut_seq=ref_mut_seq,
                                                                               gene_table=genetable,
                                                                               countinfo=countinfo, 
                                                                               Idx=idx,
                                                                               kmer=arg.kmer)
        R['output_metadata_list'], \
        R['output_kmer_lists'] = get_and_write_peptide_and_kmer(gene=gene,
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
                                                                output_silence=arg.output_silence, 
                                                                outbase=outbase)

        R['time'] = timeit.default_timer() - start_time
        R['memory'] = print_memory_diags(disable_print=True)
        R['outbase'] = outbase
        results.append(R)

    return results



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
    graph_data = graph_data[400:1400]
    remove_annot = arg.remove_annot

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
        global gene_id_list
        global sample
        global dict_kmer_foregr
        global dict_kmer_back
        global dict_pept_forgrd
        global dict_pept_backgrd
        global filepointer


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

        junction_meta_file_path = os.path.join(output_path, mutation.mode + '_metadata.tsv.gz.pq')
        gzip_tag = ''
        if arg.compressed:
            gzip_tag = '.gz'
        junction_peptide_file_path = os.path.join(output_path, mutation.mode + '_peptides.fa.pq'+gzip_tag)
        background_peptide_file_path = os.path.join(output_path, mutation.mode + '_back_peptides.fa.pq'+gzip_tag)
        junction_kmer_file_path = os.path.join(output_path, mutation.mode + '_junction_kmer.pq'+gzip_tag)
        background_kmer_file_path = os.path.join(output_path, mutation.mode + '_back_kmer.pq'+gzip_tag)
        gene_expr_file_path = os.path.join(output_path, 'gene_expression_detail.tsv'+gzip_tag)
        gene_expr_fp = gz_and_normal_open(gene_expr_file_path, 'w')
        filepointer = initialize_parquet(junction_peptide_file_path, junction_meta_file_path, background_peptide_file_path,
                    junction_kmer_file_path, background_kmer_file_path)


        dict_kmer_foregr = {}
        dict_kmer_back = {} 
        dict_pept_forgrd = {}
        dict_pept_backgrd = {}
        # go over each gene in splicegraph
        gene_id_list = list(range(0,num))
        if arg.parallel > 1:
            cnt = 0
            logging.info('Parallel: {} Threads'.format(arg.parallel)) #TODO Remove

            def process_result(gene_results):
                global cnt
                global gene_name_expr_distr
                global expr_distr
                global gene_id_list
                global sample
                global dict_kmer_foregr
                global dict_kmer_back
                global dict_pept_forgrd
                global dict_pept_backgrd
                global filepointer
                pool_genes = defaultdict(list, {})
                for gene_result in gene_results:
                    cnt += 1
                    logging.info('> Finished processing Gene {} ({}/{})'.format(gene_result['gene_name'], cnt, len(gene_id_list)))
                    if gene_result['processed']:
                        gene_name_expr_distr.append((gene_result['gene_name'], gene_result['total_expr']))
                        expr_distr.append(gene_result['total_expr'])
                        for result_type in ['background_peptide_list', 'background_kmer_lists', 'output_metadata_list', 'output_kmer_lists' ]:
                            pool_genes[result_type].extend(gene_result[result_type])
                s1 = timeit.default_timer()
                logging.debug('start writing results')
                dict_pept_forgrd, dict_pept_backgrd, dict_kmer_foregr, dict_kmer_back = write_gene_result(pool_genes, dict_pept_forgrd, dict_pept_backgrd, dict_kmer_foregr, dict_kmer_back, filepointer, remove_annot, gene_result['outbase'])
                logging.info('writing results took {} seconds'.format(timeit.default_timer() - s1))
                logging.info(">{}: {}/{} processed, time cost: {}, memory cost:{} GB ".format(sample, gene_result['gene_idx'] + 1, len(gene_id_list), gene_result['time'], gene_result['memory']))
                del gene_results

            pool = mp.Pool(processes=arg.parallel, initializer=lambda: sig.signal(sig.SIGINT, sig.SIG_IGN))
            #for gene_idx in gene_id_list:
            batch_size = 10
            for i in range(0, len(gene_id_list), batch_size):
                gene_idx = gene_id_list[i:min(i + batch_size, len(gene_id_list))]
                outbase = os.path.join(output_path, 'tmp_out_%i' % i)
                pathlib.Path(outbase).mkdir(exist_ok= True, parents= True)
                _ = pool.apply_async(process_gene_batch, args=(sample, graph_data[gene_idx], graph_info[gene_idx], gene_idx, len(gene_id_list), mutation, junction_dict, countinfo, genetable, arg, outbase,), callback=process_result)

            pool.close()
            pool.join()
            logging.debug('start collecting results')
            collect_results(filepointer.background_peptide_fp, output_path, logging)
            collect_results(filepointer.junction_peptide_fp,output_path,logging)
            collect_results(filepointer.junction_meta_fp,output_path,logging)
            collect_results(filepointer.junction_kmer_fp,output_path,logging)
            collect_results(filepointer.background_kmer_fp,output_path,logging)
            #remove_folder_list(os.path.join(output_path, 'tmp_out_')) #TODO add back after development 

        else:
            logging.info('Not Parallel')
            for gene_idx in gene_id_list:
                outbase = os.path.join(output_path, 'tmp_out')
                gene_result = process_gene_batch(sample, [graph_data[gene_idx]], [graph_info[gene_idx]], [gene_idx], len(gene_id_list), mutation, junction_dict, countinfo, genetable, arg, outbase)[0]
                if gene_result['processed']:
                    gene_name_expr_distr.append((gene_result['gene_name'], gene_result['total_expr']))
                    expr_distr.append(gene_result['total_expr'])
                    dict_pept_forgrd, dict_pept_backgrd, dict_kmer_foregr, dict_kmer_back = write_gene_result(gene_result, dict_pept_forgrd, dict_pept_backgrd, dict_kmer_foregr, dict_kmer_back,  filepointer, remove_annot, outbase)
                    time_list.append(gene_result['time'])
                    memory_list.append(gene_result['memory'])

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
            logging.info(">>>> No gene is processed during this run.")
        expr_distr_dict[sample] = expr_distr
        write_gene_expr(gene_expr_fp,gene_name_expr_distr)
        create_libsize(expr_distr_dict,output_libszie_fp)


        if arg.compressed:
            compression = 'gzip'
        else:
            compression = None


        if remove_annot:
            save_forgrd_pep_dict(dict_pept_forgrd, filepointer, compression)
            del dict_pept_forgrd

            logging.info(">>>> Foreground kmer TOTAL before final filtering {}".format(len(dict_kmer_foregr)))
            dict_kmer_foregr = filter_onkey_dict(dict_kmer_foregr, dict_kmer_back)
            logging.info(">>>> Foreground kmer TOTAL after final filtering {}".format(len(dict_kmer_foregr)))
            if len(dict_kmer_foregr) == 0 :
                logging.info(">>>> No foreground kmers are left after filtering on annotation")

            save_forgrd_kmer_dict(dict_kmer_foregr, filepointer, compression)
            del dict_kmer_foregr

            save_backgrd_kmer_dict(dict_kmer_back, filepointer, compression)
            del dict_kmer_back


        filepointer.junction_peptide_fp['filepointer'].close()
        filepointer.junction_meta_fp['filepointer'].close()
        filepointer.background_peptide_fp['filepointer'].close()
        filepointer.junction_kmer_fp['filepointer'].close()
        filepointer.background_kmer_fp['filepointer'].close()
    logging.info(">>>>>>>>> Build: Finish traversing splicegraph in mutation mode {}.\n".format(mutation.mode))

