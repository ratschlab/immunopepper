from collections import namedtuple
import glob
import gzip
import numpy as np
import os
import logging
import pyarrow as pa
import pyarrow.parquet as pq
from uuid import uuid4
import shutil
import sys
import timeit

from numpy.lib.recfunctions import structured_to_unstructured
from immunopepper.namedtuples import Filepointer

# --- intermediate fix to load pickle files stored under previous version

import spladder.classes.gene as gene
import spladder.classes.splicegraph as splicegraph
import spladder.classes.segmentgraph as segmentgraph
sys.modules['modules.classes.gene'] = gene
sys.modules['modules.classes.splicegraph'] = splicegraph
sys.modules['modules.classes.segmentgraph'] = segmentgraph


# --- end fix


def open_gz_or_normal(file_path: str, mode: str = 'r'):
    if file_path.endswith('.gz'):
        mode += 't'
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


def _convert_list_to_str(_list, sep='/'):
    remove_none_list = filter(lambda x: x is not None, _list)
    return sep.join([str(_item) for _item in remove_none_list])


def namedtuple_to_str(named_tuple: namedtuple, field_list: list[str] = None, sep: str = '\t'):
    if field_list is None:
        field_list = named_tuple._fields
    line = ''
    for field in field_list:
        if field == 'new_line':
            line = line.strip() + '\n'
            continue
        item = getattr(named_tuple, field)
        if isinstance(item, (list, tuple)):
            line += _convert_list_to_str(item, sep = ';') + sep
        else:
            line += str(item) + sep
    return line[:-1]  # remove the last separator


def list_to_tuple(list_or_tuple):
    return tuple(list_or_tuple) if isinstance(list_or_tuple, list) else list_or_tuple


def decode_utf8(s):
    return s.decode('utf-8') if hasattr(s, 'decode') else s


def get_save_path(file_info: dict, out_dir: str = None, kmer_length: int = None, create_partitions: bool = False):
    """ Parse the path stored by the filepointer depending on:
    - the nesting level of the filepointer
    - the use of a batch directory
    - the creation of parquet partitions

    :param file_info: filepointer namedtuple containing either a single path or a dictionary of paths
    :param out_dir: str. any base directory used for the path. Used for creating batches base directories
    :param kmer_length: str. the kmer length wished
    :param create_partitions: bool. whether to add an additional path depth and save uniquely identified files
    :return:

    """

    # Access the right part of the filepointer
    if kmer_length:
        path = file_info['path'][kmer_length]
    else:
        path = file_info['path']

    # Add a batch subdirectory
    if out_dir:
        path = os.path.join(out_dir, os.path.basename(path))

    # Add a parquet partition
    if create_partitions:
        path = os.path.join(path, f'part-{str(uuid4())}.pq')

    return path


def save_bg_kmer_set(data: set[str], filepointer: Filepointer, kmer_length: int, compression: str = None,
                     out_dir: str = None, verbose: bool = False):
    """ Writes a set of strings (protein kmers) to a parquet file """
    if data:
        dt = 'S' + str(kmer_length)
        data = np.fromiter(data, dt, -1) #set to array
        data = np.expand_dims(data, 0)
        path = get_save_path(filepointer.background_kmer_fp, out_dir, kmer_length)
        save_pd_toparquet(path, data, filepointer.background_kmer_fp['columns'], compression, verbose)


def save_fg_kmer_dict(data: dict, filepointer: Filepointer, kmer_length: int, compression: str = None,
                      out_dir: str = None, verbose: bool = False):
    """
    Writes kmer (of peptides) information to a parquet file.
    File structure is: ['kmer', 'id', 'segmentExpr', 'isCrossJunction', 'junctionExpr']
    """
    if data:
        data = set(namedtuple_to_str(item, sep='\t') for item in data)
        data = np.array([line.split('\t') for line in data]).T #set to array
        path = get_save_path(filepointer.junction_kmer_fp, out_dir, kmer_length)
        save_pd_toparquet(path, data, filepointer.junction_kmer_fp['columns'], compression, verbose)


def save_bg_peptide_set(data, filepointer: Filepointer, compression: str = None, out_dir: str = None,
                         verbose: bool = False):
    """
    Writes peptide information (dict of peptide to output ids) as a fasta file in parquet format.
    The keys in the fasta file are the concatenated output ids, and the value is the protein itself, e.g.::
            >ENSMUST00000027035.9
            MSSPDAGYASDDQSQPRSAQPAVMAGLGPCPWAESLSPLGDVKVKGEVVASSGAPAGTSGRAKAESRIRRPMNAF
    """
    fasta_prefix = lambda t: '>' + t
    make_fasta_ids = np.vectorize(fasta_prefix)

    if data:
        data = np.array([line.split('\t') for line in data]).T #set to array #TODO vectorize slit
        data = np.array([make_fasta_ids(data[0]), data[1]]) # output_ids,  peptides
        data = data.flatten(order = 'F')
        data = np.expand_dims(data, axis=0)
        path = get_save_path(filepointer.background_peptide_fp, out_dir)
        save_pd_toparquet(path, data, filepointer.background_peptide_fp['columns'], compression, verbose)


def save_fg_peptide_set(data: set, filepointer: Filepointer, compression: str = None, out_dir: str = None,
                         save_fasta: bool = False, verbose: bool = False):
    """
    Save foreground peptide data in a parquet file
    :param data: set containing the lines to be written sep is '\t'
    """
    fasta_prefix = lambda t: '>' + t
    make_fasta_ids = np.vectorize(fasta_prefix)

    if data:
        data = np.array([line.split('\t') for line in data]).T  # set to array
        if save_fasta:
            path_fa = get_save_path(filepointer.junction_peptide_fp, out_dir)
            fasta = np.array([make_fasta_ids(data[1]), data[0]]) #id, peptide
            fasta = fasta.flatten(order='F')
            fasta = np.expand_dims(fasta, axis=0)
            save_pd_toparquet(path_fa, fasta, filepointer.junction_peptide_fp['columns'], compression, verbose)
            del fasta

        path = get_save_path(filepointer.junction_meta_fp, out_dir)
        save_pd_toparquet(path, data, filepointer.junction_meta_fp['columns'], compression, verbose)  # Test keep peptide name in the metadata file


def save_gene_expr_distr(data: list[list], input_samples: list[str], process_sample: str, filepointer: Filepointer,
                         out_dir: str, compression: str, verbose: bool):
    """
    Save gene expression data to a parquet file that looks like this::
    +-----------------------+---------------+
    | gene                  |   ENCSR000BZG |
    |-----------------------+---------------|
    | ENSMUSG00000025902.13 |        366155 |
    | ENSMUSG00000025903.14 |        839804 |
    +-----------------------+---------------+
    :param data: list of (gene_name, gene_expression) elements
    """
    if data:
        data = np.array(data).T
        path = get_save_path(filepointer.gene_expr_fp, out_dir)
        if process_sample == 'cohort':
            data_columns = filepointer.gene_expr_fp['columns'] + list(input_samples)
        else:
            data_columns = filepointer.gene_expr_fp['columns'] + [process_sample]
        save_pd_toparquet(path, data, data_columns, compression, verbose)


def save_kmer_matrix(data_edge, data_segm, kmer_len, graph_samples, filepointer: Filepointer, compression: str = None,
                     out_dir: str = None, verbose: bool = False):
    '''
    Saves matrices of [kmers] x [metadata, samples] containing either segment expression or junction expression
    :param data_edge: set of tuples: (kmer, is_junction, junction_is_annotated, reading_frame_is_annotated, edge_expression_sample_i ..)
    :param data_segm: set of tuples: (kmer, is_junction, junction_is_annotated, reading_frame_is_annotated, segment_expression_sample_i ..)
    :param kmer_len: int amino acid length of kmer
    :param graph_samples: list sample names
    :param filepointer: class object which contains the saving path,
                        columns names and writer object for each of the files to save
    :param compression: str parquet compression format
    :param out_dir: str output directory path
    :param verbose: int verbose parameter
    '''
    logging.info("Start get file pointer")
    segm_path = get_save_path(filepointer.kmer_segm_expr_fp, out_dir, create_partitions=True)
    edge_path = get_save_path(filepointer.kmer_edge_expr_fp, out_dir, create_partitions=True)
    logging.info("Stop get file pointer")
    logging.info(f'Number of samples is {len(graph_samples)}')
    logging.info("Start creating data type array")
    dt = [('kmer', 'S' + str(kmer_len)), ('isCrossJunction', 'bool'), ('junctionAnnotated', 'bool'),
          ('readFrameAnnotated', 'bool') ] + [ (f'{sample}', 'float') for sample in graph_samples]
    logging.info("Stop creating datatype array")
    if data_edge:
        logging.info("***Start Save EDGE***")
        logging.info(f'Number of kmers to save (x patients) is {len(data_edge)}')
        logging.info("...Start From Iter")
        data_edge = np.fromiter(data_edge, dt, -1) #set to array
        logging.info("...Stop From Iter")
        logging.info("...Start structured to unstructured")
        data_edge = structured_to_unstructured(data_edge)
        logging.info("...Stop structured to unstructured")
        logging.info("...Start transpose")
        data_edge = data_edge.T
        logging.info("...Stop transpose")
        data_edge_columns = filepointer.kmer_segm_expr_fp['columns'] + graph_samples
        filepointer.kmer_edge_expr_fp['pqwriter'] = save_pd_toparquet(edge_path, data_edge, data_edge_columns,
                                                                      compression=compression, verbose=verbose,
                                                                      pqwriter=filepointer.kmer_edge_expr_fp[
                                                                          'pqwriter'], writer_close=True)
        logging.info("***Stop Save EDGE***")
    if data_segm:
        logging.info("****Start save SEGM***")
        logging.info(f'Number of kmers to save (x patients) is {len(data_segm)}')
        logging.info("...Start From Iter")
        data_segm = np.fromiter(data_segm, dt, -1) #set to array
        logging.info("...Stop From Iter")
        logging.info("...Start structured to unstructured")
        data_segm = structured_to_unstructured(data_segm)
        logging.info("...Stop structured to unstructured")
        logging.info("...Start transpose")
        data_segm = data_segm.T
        logging.info("...Stop transpose")
        data_segm_columns = filepointer.kmer_segm_expr_fp['columns'] + graph_samples
        filepointer.kmer_segm_expr_fp['pqwriter'] = save_pd_toparquet(segm_path, data_segm, data_segm_columns,
                                                                      compression=compression, verbose=verbose,
                                                                      pqwriter=filepointer.kmer_segm_expr_fp[
                                                                          'pqwriter'], writer_close=True)
        logging.info("****Stop save SEGM****")




def initialize_fp(output_path: str, mutation_mode: str,
                  kmer_list: list[int], output_fasta: bool, cross_graph_expr: bool):
    """"
    Initializes a Filepointer object for the given parameters.
    :param output_path: base directory where data is written
    :param mutation_mode: what type of mutation mode we operate in, such as 'ref', 'somatic', 'germline', etc.
    :kmer_list: list of kmer sizes to operate on
    :output_fasta: if true, create a fasta file for foreground peptides
    :cross_graph_expr: if true, write gene expression matrices with counts from the full graph
    :return: Filepointer instance containing output information for all information that will be written to a file
    :rtype: Filepointer
    """

    # --- Paths
    annot_peptide_file_path = os.path.join(output_path, mutation_mode + '_annot_peptides.fa.pq')
    annot_kmer_file_path = os.path.join(output_path, mutation_mode + '_annot_kmer.pq')
    gene_expr_file_path = os.path.join(output_path, 'gene_expression_detail.pq')
    junction_meta_file_path = os.path.join(output_path, mutation_mode + '_sample_peptides_meta.pq')
    junction_peptide_file_path = os.path.join(output_path, mutation_mode + '_sample_peptides.fa.pq')
    junction_kmer_file_path = os.path.join(output_path, mutation_mode + '_sample_kmer.pq')
    graph_kmer_segment_expr_path = os.path.join(output_path, mutation_mode + '_graph_kmer_SegmExpr')
    graph_kmer_junction_expr_path = os.path.join(output_path, mutation_mode + '_graph_kmer_JuncExpr')

    # --- Fields
    cols_annot_pep_file = ['fasta']
    cols_annot_kmer_file = ['kmer']
    cols_gene_expr_file = ['gene']
    cols_metadata_file = ['peptide', 'id', 'readFrame', 'readFrameAnnotated', 'geneName', 'geneChr', 'geneStrand',
                                'mutationMode',
                                'junctionAnnotated', 'hasStopCodon', 'isInJunctionList',
                                'isIsolated', 'variantComb', 'variantSegExpr', 'modifiedExonsCoord',
                                'originalExonsCoord', 'vertexIdx', 'junctionExpr', 'segmentExpr',
                                'kmerType']
    cols_pep_file = ['fasta']
    cols_kmer_single_sample_file = ['kmer', 'id', 'segmentExpr', 'isCrossJunction', 'junctionExpr',
                               'junctionAnnotated', 'readFrameAnnotated']
    cols_metaOnly_kmer_cross_sample_file = ['kmer', 'isCrossJunction', 'junctionAnnotated', 'readFrameAnnotated']

    # --- Grouping dict
    if output_fasta:  # Foreground peptide fasta - optional
        peptide_fp = _output_info(junction_peptide_file_path, cols_pep_file)
    else:
        peptide_fp = None
    if cross_graph_expr:  # Expression matrices with counts from full graph - optional
        junction_kmer_fp = None
        kmer_segm_expr_fp = _output_info(graph_kmer_segment_expr_path, cols_metaOnly_kmer_cross_sample_file, pq_writer=True)
        kmer_edge_expr_fp = _output_info(graph_kmer_junction_expr_path, cols_metaOnly_kmer_cross_sample_file, pq_writer=True)
    else:  # Expression kmer information from single sample
        junction_kmer_fp = _output_info(junction_kmer_file_path, cols_kmer_single_sample_file, kmer_list)
        kmer_segm_expr_fp = None
        kmer_edge_expr_fp = None
    metadata_fp = _output_info(junction_meta_file_path, cols_metadata_file)
    annot_fp = _output_info(annot_peptide_file_path, cols_annot_pep_file)
    annot_kmer_fp = _output_info(annot_kmer_file_path, cols_annot_kmer_file, kmer_list)
    gene_expr_fp = _output_info(gene_expr_file_path, cols_gene_expr_file)

    # --- Filepointer object
    filepointer = Filepointer(peptide_fp, metadata_fp, annot_fp, junction_kmer_fp, annot_kmer_fp,
                              gene_expr_fp, kmer_segm_expr_fp, kmer_edge_expr_fp)
    return filepointer


def _output_info(path, file_columns, kmer_lengths=None, pq_writer=False):
    """ Builds a path for each kmer length """
    if kmer_lengths:  # variable number of kmer files
        path_dict = {}
        for kmer_length in kmer_lengths:
            path_dict[kmer_length] = path.replace('kmer.pq', f'{kmer_length}mer.pq')
        path = path_dict
    file_info = {'path': path, 'columns': file_columns}
    if pq_writer:
        file_info['pqwriter'] = None
    return file_info


def save_pd_toparquet(path, array, columns, compression=None, verbose=False, pqwriter=None, writer_close=True):
    """
    Saves a pandas data frame in parquet format.
    """
    s1 = timeit.default_timer()
    logging.info("......Start convert pyarrow")
    table = pa.Table.from_arrays(array, names = columns)
    logging.info("......Stop convert pyarrow")
    if pqwriter is None:
        logging.info("......Start writer")
        pqwriter = pq.ParquetWriter(path, table.schema, compression=compression)
    pqwriter.write_table(table)
    logging.info("......Stop writer")
    if writer_close:
        logging.info("......Start close writer")
        pqwriter.close()
        logging.info("......Stop close writer")
        pqwriter = None

    if verbose:
        file_name = os.path.basename(path)
        tot_shape = array.shape[1]
        logging.info(f'Saved {file_name} with {tot_shape} lines in {round(timeit.default_timer() - s1, 4)}s')
    return pqwriter


def collect_results(filepointer_item, out_dir, compression, mutation_mode, kmer_list=None, parquet_partitions=False):
    """
    Merges results written by each parallel process into a single rechunked parquet file.
    """
    if filepointer_item is not None:
        # get paths for several kmer lengths
        if kmer_list:
            file_to_collect = filepointer_item['path'].values()
        else:
            file_to_collect = [filepointer_item['path']]

        # merge the partitions
        for file_path in file_to_collect:
            file_name = os.path.basename(file_path)
            # adjust name for parquet partitions
            if parquet_partitions:
                tmp_file_list = glob.glob(os.path.join(out_dir, f'tmp_out_{mutation_mode}_batch_[0-9]*',
                                                       file_name, '*part*'))
            else:
                tmp_file_list = glob.glob(os.path.join(out_dir, f'tmp_out_{mutation_mode}_batch_[0-9]*',
                                                       file_name))

            try:
                dset = pq.ParquetDataset(tmp_file_list)
                dset_pandas = dset.read().to_pandas() #TODO check efficiency
                save_pd_toparquet(file_path, dset_pandas, compression, verbose=1)
            except:
                logging.error(f'Unable to read one of files for {file_path} collection')
                sys.exit(1)


def remove_folder_list(base_path):
    """ Removes folders starting with base_path. Used for cleaning up temporary files in multiprocessing mode. """
    folder_list = glob.glob(base_path + '*')
    for dir_path in folder_list:
        shutil.rmtree(dir_path, ignore_errors=True)


def read_pq_with_dict(file_path, columns):
    return pa.parquet.read_table(file_path, read_dictionary=columns)


def read_pq_with_pd(file_path):
    return pa.parquet.read_table(file_path).to_pandas()
