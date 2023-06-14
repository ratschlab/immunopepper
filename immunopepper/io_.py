from collections import namedtuple
import csv
import glob
import gzip
import numpy as np
import os
import logging
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from uuid import uuid4
import shutil
import sys
from time import sleep
import timeit

from immunopepper.namedtuples import Filepointer

# --- intermediate fix to load pickle files stored under previous version

import spladder.classes.gene as gene
import spladder.classes.splicegraph as splicegraph
import spladder.classes.segmentgraph as segmentgraph
sys.modules['modules.classes.gene'] = gene
sys.modules['modules.classes.splicegraph'] = splicegraph
sys.modules['modules.classes.segmentgraph'] = segmentgraph


# --- end fix

fasta_prefix = lambda t: '>' + t
make_fasta_ids = np.vectorize(fasta_prefix)


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


def get_save_path(file_info: dict, out_dir: str = None,  create_partitions: bool = False):
    """ Parse the path stored by the filepointer depending on:
    - the nesting level of the filepointer
    - the use of a batch directory
    - the creation of parquet partitions

    :param file_info: filepointer namedtuple containing either a single path or a dictionary of paths
    :param out_dir: str. any base directory used for the path. Used for creating batches base directories
    :param create_partitions: bool. whether to add an additional path depth and save uniquely identified files
    :return:

    """

    # Access the right part of the filepointer
    path = file_info['path']

    # Add a batch subdirectory
    if out_dir:
        path = os.path.join(out_dir, os.path.basename(path))

    # Add a parquet partition
    if create_partitions:
        path = os.path.join(path, f'part-{str(uuid4())}.gz')

    return path


def save_bg_kmer_set(data: set[str], filepointer: Filepointer,
                     out_dir: str = None, verbose: bool = False):
    """ Writes a set of strings (protein kmers)""" #TODO update docs
    if data:
        path = get_save_path(filepointer.background_kmer_fp, out_dir)
        save_to_gzip(path, data, filepointer.background_kmer_fp['columns'], verbose=verbose)


def save_fg_kmer_dict(data: dict, filepointer: Filepointer, kmer_length: int,
                      out_dir: str = None, verbose: bool = False):
    """
    Writes kmer (of peptides) information.
    File structure is: ['kmer', 'id', 'segmentExpr', 'isCrossJunction', 'junctionExpr']
    """
    if data:
        data = set(namedtuple_to_str(item, sep='\t') for item in data)
        path = get_save_path(filepointer.junction_kmer_fp, out_dir, kmer_length)
        save_to_gzip(path, data, filepointer.junction_kmer_fp['columns'], verbose)


def save_bg_peptide_set(data, filepointer: Filepointer, out_dir: str = None,
                         verbose: bool = False):
    """
    Writes peptide information (dict of peptide to output ids) as a fasta file.
    The keys in the fasta file are the concatenated output ids, and the value is the protein itself, e.g.::
            >ENSMUST00000027035.9
            MSSPDAGYASDDQSQPRSAQPAVMAGLGPCPWAESLSPLGDVKVKGEVVASSGAPAGTSGRAKAESRIRRPMNAF
    """

    if data:
        data = np.array([line.split('\t') for line in data]).T #set to array #TODO vectorize slit
        data = np.array([make_fasta_ids(data[0]), data[1]]) # output_ids,  peptides
        data = data.flatten(order = 'F')
        data = np.expand_dims(data, axis=0).T
        path = get_save_path(filepointer.background_peptide_fp, out_dir)
        save_to_gzip(path, data, filepointer.background_peptide_fp['columns'], verbose, is_2d=True)


def save_fg_peptide_set(data: set, filepointer: Filepointer, out_dir: str = None,
                         save_fasta: bool = False, verbose: bool = False):
    """
    Save foreground peptide data.
    :param data: set containing the lines to be written sep is '\t'
    """

    if data:
        data = np.array([line.split('\t') for line in data]).T  # set to array
        if save_fasta:
            path_fa = get_save_path(filepointer.junction_peptide_fp, out_dir)
            fasta = np.array([make_fasta_ids(data[1]), data[0]]) #id, peptide
            fasta = fasta.flatten(order='F')
            fasta = np.expand_dims(fasta, axis=0).T
            save_to_gzip(path_fa, fasta, filepointer.junction_peptide_fp['columns'], verbose, is_2d=True)
            del fasta

        path = get_save_path(filepointer.junction_meta_fp, out_dir)
        save_to_gzip(path, data.T, filepointer.junction_meta_fp['columns'], verbose, is_2d=True)  # Test keep peptide name in the metadata file


def save_gene_expr_distr(data: list[list], input_samples: list[str], process_sample: str, filepointer: Filepointer,
                         out_dir: str, verbose: bool):
    """
    Save gene expression data to a file that looks like this::
    +-----------------------+---------------+
    | gene                  |   ENCSR000BZG |
    |-----------------------+---------------|
    | ENSMUSG00000025902.13 |        366155 |
    | ENSMUSG00000025903.14 |        839804 |
    +-----------------------+---------------+
    :param data: list of (gene_name, gene_expression) elements
    """
    if data:
        path = get_save_path(filepointer.gene_expr_fp, out_dir)
        if process_sample == 'cohort':
            data_columns = filepointer.gene_expr_fp['columns'] + list(input_samples)
        else:
            data_columns = filepointer.gene_expr_fp['columns'] + [process_sample]
        save_to_gzip(path, data, data_columns, verbose=verbose, is_2d=True)


def save_kmer_matrix(data_edge, data_segm, graph_samples, filepointer: Filepointer,
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
    segm_path = get_save_path(filepointer.kmer_segm_expr_fp, out_dir, create_partitions=True)
    edge_path = get_save_path(filepointer.kmer_edge_expr_fp, out_dir, create_partitions=True)
    if data_edge:
        data_edge_columns = filepointer.kmer_segm_expr_fp['columns'] + graph_samples
        filepointer.kmer_edge_expr_fp['pqwriter'] = save_to_gzip(edge_path, data_edge, data_edge_columns,
                                                                 verbose=verbose,
                                                                 filepointer=filepointer.kmer_edge_expr_fp['pqwriter'],
                                                                 writer_close=True, is_2d=True)
    if data_segm:
        data_segm_columns = filepointer.kmer_segm_expr_fp['columns'] + graph_samples
        filepointer.kmer_segm_expr_fp['pqwriter'] = save_to_gzip(segm_path, data_segm, data_segm_columns,
                                                                 verbose=verbose,
                                                                 filepointer=filepointer.kmer_segm_expr_fp[ 'pqwriter'],
                                                                 writer_close=True, is_2d=True)




def initialize_fp(output_path: str, mutation_mode: str, output_fasta: bool):
    """"
    Initializes a Filepointer object for the given parameters.
    :param output_path: base directory where data is written
    :param mutation_mode: what type of mutation mode we operate in, such as 'ref', 'somatic', 'germline', etc.
    :output_fasta: if true, create a fasta file for foreground peptides
    :return: Filepointer instance containing output information for all information that will be written to a file
    :rtype: Filepointer
    """

    # --- Paths
    annot_peptide_file_path = os.path.join(output_path, mutation_mode + '_annot_peptides.fa.gz')
    annot_kmer_file_path = os.path.join(output_path, mutation_mode + '_annot_kmer.gz')
    gene_expr_file_path = os.path.join(output_path, 'gene_expression_detail.gz')
    junction_meta_file_path = os.path.join(output_path, mutation_mode + '_sample_peptides_meta.gz')
    junction_peptide_file_path = os.path.join(output_path, mutation_mode + '_sample_peptides.fa.gz')
    graph_kmer_segment_expr_path = os.path.join(output_path, mutation_mode + '_graph_kmer_SegmExpr')
    graph_kmer_junction_expr_path = os.path.join(output_path, mutation_mode + '_graph_kmer_JuncExpr')

    # --- Fields
    cols_annot_pep_file = ['fasta']
    cols_annot_kmer_file = ['kmer']
    cols_gene_expr_file = ['gene']
    cols_metadata_file = ['peptide', 'id', 'readFrame', 'readFrameAnnotated', 'geneName', 'geneChr', 'geneStrand',
                          'mutationMode', 'hasStopCodon', 'isInJunctionList', 'isIsolated', 'variantComb',
                          'variantSegExpr', 'modifiedExonsCoord', 'originalExonsCoord', 'vertexIdx', 'kmerType']
    cols_pep_file = ['fasta']
    metacols_kmers_expr_file = ['kmer', 'coord', 'isCrossJunction', 'junctionAnnotated', 'readFrameAnnotated']

    # --- Grouping dict
    if output_fasta:  # Foreground peptide fasta - optional
        peptide_fp = _output_info(junction_peptide_file_path, cols_pep_file)
    else:
        peptide_fp = None

    kmer_segm_expr_fp = _output_info(graph_kmer_segment_expr_path, metacols_kmers_expr_file, pq_writer=True)
    kmer_edge_expr_fp = _output_info(graph_kmer_junction_expr_path, metacols_kmers_expr_file, pq_writer=True)

    metadata_fp = _output_info(junction_meta_file_path, cols_metadata_file)
    annot_fp = _output_info(annot_peptide_file_path, cols_annot_pep_file)
    annot_kmer_fp = _output_info(annot_kmer_file_path, cols_annot_kmer_file)
    gene_expr_fp = _output_info(gene_expr_file_path, cols_gene_expr_file)

    # --- Filepointer object
    filepointer = Filepointer(peptide_fp, metadata_fp, annot_fp, annot_kmer_fp,
                              gene_expr_fp, kmer_segm_expr_fp, kmer_edge_expr_fp)
    return filepointer


def _output_info(path, file_columns, pq_writer=False):
    """ Creates a dictionary with keys:
    path: the path of the file
    columns: the columns of the file
    pqwriter: the file handle """
    file_info = {'path': path, 'columns': file_columns}
    if pq_writer:
        file_info['pqwriter'] = None
    return file_info


def disk_writer(path, data_iterable, columns, filepointer=None, writer_close=True, is_2d=False):
    """Write gzip to disk"""
    delim = '\t'
    linet = '\n'
    # Open
    if filepointer is None:
        filepointer = gzip.open(path, 'wt')
        filepointer.write(delim.join(columns) + linet)

    # Write
    if is_2d:
        for idx, line in enumerate(data_iterable):
            filepointer.write(delim.join([str(x) for x in line]) + linet)
    else:
        for line in data_iterable:
            filepointer.write(line + '\n')

    # Close
    if writer_close:
        filepointer.close()
        filepointer = None


def save_to_gzip(path, data_iterable, columns, verbose=False, filepointer=None, writer_close=True, is_2d=False):
    """
    Saves a  data frame in gzip format.
    """
    s1 = timeit.default_timer()

    error_encountered = 1
    sleep_times = [60, 10, 1]  # Number of seconds for each saving attempt
    while sleep_times and error_encountered:
        try:
            disk_writer(path, data_iterable, columns, filepointer, writer_close, is_2d)
            error_encountered = 0
        except OSError:
            sleep_time = sleep_times.pop()
            logging.info(f'Issue saving {path} to file system: sleeping {sleep_time} sec. and retry')
            sleep(sleep_time)
            error_encountered = 1

    if error_encountered:
        raise OSError('Issue with saving device')

    if verbose:
        file_name = os.path.basename(path)
        tot_shape = len(data_iterable)
        logging.info(f'Saved {file_name} with {tot_shape} lines in {round(timeit.default_timer() - s1, 4)}s')
    return filepointer


def collect_results(filepointer_item, out_dir, mutation_mode, partitions=False):
    """
    Merges results written by each parallel process.
    """
    if filepointer_item is not None:
        file_path = filepointer_item['path']

    # merge the partitions
    file_name = os.path.basename(file_path)
    # adjust name for partitions
    if partitions:
        tmp_file_list = glob.glob(os.path.join(out_dir, f'tmp_out_{mutation_mode}_batch_[0-9]*', file_name, '*part*'))
    else:
        tmp_file_list = glob.glob(os.path.join(out_dir, f'tmp_out_{mutation_mode}_batch_[0-9]*', file_name))

    try:
        df = pd.concat((pd.read_csv(f, sep='\t', compression='gzip') for f in tmp_file_list), ignore_index=True)
    except:
        logging.error(f'Unable to read one of files for {file_path} collection')
        sys.exit(1)
    return df


def remove_folder_list(base_path):
    """ Removes folders starting with base_path. Used for cleaning up temporary files in multiprocessing mode. """
    folder_list = glob.glob(base_path + '*')
    for dir_path in folder_list:
        shutil.rmtree(dir_path, ignore_errors=True)


def read_pq_with_dict(file_path, columns):
    return pa.parquet.read_table(file_path, read_dictionary=columns)


def read_pq_with_pd(file_path):
    return pa.parquet.read_table(file_path).to_pandas()
