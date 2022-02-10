from collections import namedtuple
import glob
import gzip
import os
import logging
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import shutil
import sys
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
            line += _convert_list_to_str(item) + sep
        else:
            line += str(item) + sep
    return line[:-1]  # remove the last separator


def list_to_tuple(list_or_tuple):
    return tuple(list_or_tuple) if isinstance(list_or_tuple, list) else list_or_tuple


def decode_utf8(s):
    return s.decode('utf-8') if hasattr(s, 'decode') else s


def get_save_path(file_info: dict, out_dir: str = None, kmer_length: int = None):
    """ Computes the path where data should be saved for the given file_info and kmer_length """
    if kmer_length:
        input_path = file_info['path'][kmer_length]
    else:
        input_path = file_info['path']
    if out_dir:
        return os.path.join(out_dir, os.path.basename(input_path))
    return input_path


def save_bg_kmer_set(data: set[str], filepointer: Filepointer, kmer_length: int, compression: str = None,
                     out_dir: str = None, verbose: bool = False):
    """ Writes a set of strings (protein kmers) to a parquet file """
    if data:
        data = pd.DataFrame(data, columns=['kmer'])
        path = get_save_path(filepointer.background_kmer_fp, out_dir, kmer_length)
        save_pd_toparquet(path, data, compression, verbose)


def save_fg_kmer_dict(data: dict, filepointer: Filepointer, kmer_length: int, compression: str = None,
                      out_dir: str = None, verbose: bool = False):
    """
    Writes kmer (of peptides) information to a parquet file.
    File structure is: ['kmer', 'id', 'segmentExpr', 'isCrossJunction', 'junctionExpr']
    """
    if data:
        data = pd.DataFrame(data.values(), index=data.keys())
        data = data.applymap(_convert_list_to_str)
        data = data.rename_axis('kmer').reset_index()
        data.columns = filepointer.junction_kmer_fp['columns']
        path = get_save_path(filepointer.junction_kmer_fp, out_dir, kmer_length)
        save_pd_toparquet(path, data, compression, verbose)


def save_bg_peptide_dict(data, filepointer: Filepointer, compression: str = None, out_dir: str = None,
                         verbose: bool = False):
    """
    Writes peptide information (dict of peptide to output ids) as a fasta file in parquet format.
    The keys in the fasta file are the concatenated output ids, and the value is the protein itself, e.g.::
            >ENSMUST00000027035.9/ENSMUST00000116652.7
            MSSPDAGYASDDQSQPRSAQPAVMAGLGPCPWAESLSPLGDVKVKGEVVASSGAPAGTSGRAKAESRIRRPMNAF
    """
    if data:
        data = pd.DataFrame(data.values(), index=data.keys()).reset_index()
        data.columns = ['index', 'outputId']
        data['outputId'] = data['outputId'].apply(lambda x: '>' + '/'.join(x))
        data = pd.concat([data['outputId'], data['index']]).sort_index(kind="mergesort").reset_index(drop=True)
        data = pd.DataFrame(data, columns=filepointer.background_peptide_fp['columns'])
        path = get_save_path(filepointer.background_peptide_fp, out_dir)
        save_pd_toparquet(path, data, compression, verbose)


def save_fg_peptide_dict(data: dict, filepointer: Filepointer, compression: str = None, out_dir: str = None,
                         save_fasta: bool = False, verbose: bool = False):
    """
    Save foreground peptide data in a parquet file
    :param data: dictionary mapping peptides to a list of feature sets
    """
    if data:
        data = pd.DataFrame(data.values(), index=data.keys()).reset_index()
        data.columns = filepointer.junction_meta_fp['columns']
        if save_fasta:
            path_fa = get_save_path(filepointer.junction_peptide_fp, out_dir)
            fasta = pd.concat([data['id'].apply(lambda x: '>' + '/'.join(x)), data['peptide']]).sort_index(
                kind="mergesort").reset_index(drop=True)
            fasta = pd.DataFrame(fasta, columns=filepointer.junction_peptide_fp['columns'])
            save_pd_toparquet(path_fa, fasta, compression, verbose)
            del fasta
        data = data.set_index('peptide')
        data = data.applymap(_convert_list_to_str)
        data = data.reset_index(drop=False)
        path = get_save_path(filepointer.junction_meta_fp, out_dir)
        save_pd_toparquet(path, data, compression, verbose)  # Test keep peptide name in the metadata file


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
        data = pd.DataFrame(data)
        path = get_save_path(filepointer.gene_expr_fp, out_dir)
        if process_sample == 'cohort':
            data.columns = filepointer.gene_expr_fp['columns'] + list(input_samples)
        else:
            data.columns = filepointer.gene_expr_fp['columns'] + [process_sample]
        save_pd_toparquet(path, data, compression, verbose)


def save_kmer_matrix(data, graph_samples, filepointer: Filepointer, compression: str = None, out_dir: str = None,
                     verbose: bool = False):
    """
    Saves the kmer segment and junction matrices per gene.
    The filepointer is kept open until all genes from the batch are written
    """
    if data[0]:
        segm_path = get_save_path(filepointer.kmer_segm_expr_fp, out_dir)
        edge_path = get_save_path(filepointer.kmer_edge_expr_fp, out_dir)
        data[1] = pd.DataFrame.from_dict(data[1], orient='index', columns=graph_samples).reset_index().rename(
            {'index': filepointer.kmer_segm_expr_fp['columns'][0]}, axis=1)
        data[2] = pd.DataFrame.from_dict(data[2], orient='index', columns=graph_samples).reset_index().rename(
            {'index': filepointer.kmer_segm_expr_fp['columns'][0]}, axis=1)

        data[1][filepointer.kmer_segm_expr_fp['columns'][1]] = data[0].values()
        data[2][filepointer.kmer_segm_expr_fp['columns'][1]] = data[0].values()
        filepointer.kmer_segm_expr_fp['pqwriter'] = save_pd_toparquet(segm_path, data[1],
                                                                      compression=compression, verbose=verbose,
                                                                      pqwriter=filepointer.kmer_segm_expr_fp[
                                                                          'pqwriter'], writer_close=False)
        filepointer.kmer_edge_expr_fp['pqwriter'] = save_pd_toparquet(edge_path, data[2],
                                                                      compression=compression, verbose=verbose,
                                                                      pqwriter=filepointer.kmer_edge_expr_fp[
                                                                          'pqwriter'], writer_close=False)


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
    background_peptide_file_path = os.path.join(output_path, mutation_mode + '_annot_peptides.fa.pq')
    background_kmer_file_path = os.path.join(output_path, mutation_mode + '_annot_kmer.pq')
    gene_expr_file_path = os.path.join(output_path, 'gene_expression_detail.pq')
    junction_meta_file_path = os.path.join(output_path, mutation_mode + '_sample_peptides_meta.pq')
    junction_peptide_file_path = os.path.join(output_path, mutation_mode + '_sample_peptides.fa.pq')
    junction_kmer_file_path = os.path.join(output_path, mutation_mode + '_sample_kmer.pq')
    graph_kmer_segment_expr_path = os.path.join(output_path, mutation_mode + '_graph_kmer_SegmExpr.pq')
    graph_kmer_junction_expr_path = os.path.join(output_path, mutation_mode + '_graph_kmer_JuncExpr.pq')

    # --- Fields
    fields_backgrd_pep_dict = ['fasta']
    fields_backgrd_kmer_dict = ['kmer']
    fields_gene_expr_file = ['gene']
    fields_meta_peptide_dict = ['peptide', 'id', 'readFrame', 'geneName', 'geneChr', 'geneStrand',
                                'mutationMode',
                                'junctionAnnotated', 'hasStopCodon', 'isInJunctionList',
                                'isIsolated', 'variantComb', 'variantSegExpr', 'modifiedExonsCoord',
                                'originalExonsCoord', 'vertexIdx',
                                'kmerType']
    fields_forgrd_pep_dict = ['fasta']
    fields_forgrd_kmer_dict = ['kmer', 'id', 'segmentExpr', 'isCrossJunction', 'junctionExpr']
    fields_kmer_expr = ['kmer', 'isCrossJunction']

    # --- Grouping dict
    if output_fasta:  # Foreground peptide fasta - optional
        peptide_fp = _output_info(junction_peptide_file_path, fields_forgrd_pep_dict)
    else:
        peptide_fp = None
    if cross_graph_expr:  # Expression matrices with counts from full graph - optional
        junction_kmer_fp = None
        kmer_segm_expr_fp = _output_info(graph_kmer_segment_expr_path, fields_kmer_expr, pq_writer=True)
        kmer_edge_expr_fp = _output_info(graph_kmer_junction_expr_path, fields_kmer_expr, pq_writer=True)
    else:  # Expression kmer information from single sample
        fields_meta_peptide_dict.insert(-1, 'junctionExpr')
        fields_meta_peptide_dict.insert(-1, 'segmentExpr')
        junction_kmer_fp = _output_info(junction_kmer_file_path, fields_forgrd_kmer_dict, kmer_list)
        kmer_segm_expr_fp = None
        kmer_edge_expr_fp = None
    meta_peptide_fp = _output_info(junction_meta_file_path, fields_meta_peptide_dict)
    background_fp = _output_info(background_peptide_file_path, fields_backgrd_pep_dict)
    background_kmer_fp = _output_info(background_kmer_file_path, fields_backgrd_kmer_dict, kmer_list)
    gene_expr_fp = _output_info(gene_expr_file_path, fields_gene_expr_file)

    # --- Filepointer object
    filepointer = Filepointer(peptide_fp, meta_peptide_fp, background_fp, junction_kmer_fp, background_kmer_fp,
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


def save_pd_toparquet(path, pd_df, compression=None, verbose=False, pqwriter=None, writer_close=True):
    """
    Saves a panda data frame in parquet format.
    """
    s1 = timeit.default_timer()
    size_init = pd_df.shape[1]
    pd_df = pd_df.loc[:, ~pd_df.columns.duplicated(keep='first')]
    size_new = pd_df.shape[1]
    if size_init != size_new:
        logging.warning("Graph/Count file contains duplicated sample names. Kept the first occurence")
    table = pa.Table.from_pandas(pd_df, preserve_index=False)

    if pqwriter is None:
        pqwriter = pq.ParquetWriter(path, table.schema, compression=compression)
    pqwriter.write_table(table)
    if writer_close:
        pqwriter.close()

    if verbose:
        file_name = os.path.basename(path)
        tot_shape = pd_df.shape[0]
        logging.info(f'Saved {file_name} with {tot_shape} lines in {round(timeit.default_timer() - s1, 4)}s')
    return pqwriter


def collect_results(filepointer_item, out_dir, compression, mutation_mode, kmer_list=None):
    """
    Merges results written by each parallel process into a single parquet file.
    """
    if filepointer_item is not None:
        s1 = timeit.default_timer()

        if kmer_list:
            file_to_collect = filepointer_item['path'].values()
        else:
            file_to_collect = [filepointer_item['path']]
        for file_path in file_to_collect:
            file_name = os.path.basename(file_path)
            tmp_file_list = glob.glob(os.path.join(out_dir, f'tmp_out_{mutation_mode}_batch_[0-9]*', file_name))
            tot_shape = 0
            for tmp_file in tmp_file_list:
                try:
                    table = pq.read_table(tmp_file)
                    if tot_shape == 0:
                        pqwriter = pq.ParquetWriter(file_path, table.schema, compression=compression)
                    pqwriter.write_table(table)
                    tot_shape += table.shape[0]
                except:
                    logging.error(f'Unable to read: {tmp_file}')
                    sys.exit(1)
            if tmp_file_list:
                pqwriter.close()
                logging.info(f'Collected {file_name} with {tot_shape} lines in {timeit.default_timer() - s1}s');


def remove_folder_list(base_path):
    """ Removes folders starting with base_path. Used for cleaning up temporary files in multiprocessing mode. """
    folder_list = glob.glob(base_path + '*')
    for dir_path in folder_list:
        shutil.rmtree(dir_path, ignore_errors=True)


def read_pq_with_dict(file_path, columns):
    return pa.parquet.read_table(file_path, read_dictionary=columns)


def read_pq_with_pd(file_path):
    return pa.parquet.read_table(file_path).to_pandas()
