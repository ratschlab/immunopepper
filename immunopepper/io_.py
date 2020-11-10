import glob
import gzip
import os
import logging
import pandas as pd
import pickle
import pyarrow as pa
import pyarrow.parquet as pq
import shutil
import sys
import timeit

from .namedtuples import Coord
from .namedtuples import Filepointer

### intermediate fix to load pickle files stored under previous version
import spladder.classes.gene as gene
import spladder.classes.splicegraph as splicegraph
import spladder.classes.segmentgraph as segmentgraph
sys.modules['modules.classes.gene'] = gene
sys.modules['modules.classes.splicegraph'] = splicegraph
sys.modules['modules.classes.segmentgraph'] = segmentgraph
### end fix


def load_pickled_graph(f):
    return pickle.load(f)


def gz_and_normal_open(file_path,mode='r'):
    if file_path.endswith('.gz'):
        mode += 't'
        file_fp = gzip.open(file_path, mode)
    else:
        file_fp = open(file_path, mode)

    return file_fp


def write_list(fp, _list):
    fp.writelines([l + '\n' for l in _list])


def _convert_list_to_str(_list, sep ='/'):
    remove_none_list = filter(lambda x:x is not None, _list)
    return sep.join([str(_item) for _item in remove_none_list])

def convert_namedtuple_to_str(_namedtuple, field_list = None, sep = '\t'):

    if field_list is None:
        field_list = _namedtuple._fields
    line = ''
    for field in field_list:
        if field == 'new_line':
            line = line.strip() + '\n'
            continue
        # if not hasattr(_namedtuple, field):
        #     logging.error('Namedtuple %s ' % str(_namedtuple) + ' does not have a field: ' + field)
        #     sys.exit(1)
        item = getattr(_namedtuple, field)
        if isinstance(item, (list, tuple)):
            line += _convert_list_to_str(item)+sep
        else:
            line += str(item)+sep
    return line[:-1] # remove the last '\t'



def list_to_tuple(input):
    if isinstance(input, list):
        return tuple(input)
    else:
        return input

def convert_to_str_Coord_namedtuple(input, sep = ';'):
    if isinstance(input, Coord):
        input = convert_namedtuple_to_str(input, sep = sep)
    return input


def write_namedtuple_list(fp, namedtuple_list, field_list):
    """ Write namedtuple_list to the given file pointer"""
    fp.writelines(convert_namedtuple_to_str(_namedtuple, field_list) + '\n' for _namedtuple in namedtuple_list)

def write_gene_expr(fp, gene_expr_tuple_list):
    header_line = 'gene\texpr\n'
    fp.write(header_line)
    gene_expr_str_list = [ gene_expr_tuple[0]+'\t'+str(gene_expr_tuple[1]) for gene_expr_tuple in gene_expr_tuple_list]
    write_list(fp,gene_expr_str_list)

def codeUTF8(s):
    return s.encode('utf-8')

def decodeUTF8(s):
    if not hasattr(s, 'decode'):
        return s
    return s.decode('utf-8')


def switch_tmp_path(file_info, outbase=None, kmer_length = None):
    if kmer_length:
        input_path = file_info['path'][kmer_length]
    else:
        input_path = file_info['path']
    '''switches between temporary output folder in parallel mode and main output folder'''
    if outbase:
        path = os.path.join(outbase, os.path.basename(input_path))
    else:
        path = input_path
    return path


def save_backgrd_kmer_set(data, filepointer, kmer_length, compression=None, outbase=None, verbose=False):
    path = switch_tmp_path(filepointer.background_kmer_fp, outbase, kmer_length)
    if data:
        data = pd.DataFrame(data, columns=['kmer'])
        save_pd_toparquet(path, data, compression, verbose)


def save_forgrd_kmer_dict(data, filepointer, kmer_length, compression=None, outbase=None, verbose=False):
    path = switch_tmp_path(filepointer.junction_kmer_fp, outbase, kmer_length)
    if data:
        data = pd.DataFrame(data.values(), index=data.keys())
        data = data.applymap(_convert_list_to_str)
        data = data.rename_axis('kmer').reset_index()
        data.columns = filepointer.junction_kmer_fp['columns']
        save_pd_toparquet(path, data, compression, verbose)


def save_backgrd_pep_dict(data, filepointer, compression = None, outbase=None, verbose=False):
    path = switch_tmp_path(filepointer.background_peptide_fp, outbase)
    if data:
        data = pd.DataFrame(data.values(), index = data.keys()).reset_index()
        data.columns = ['index', 'output_id']
        data['output_id'] = data['output_id'].apply(lambda x: '>' + '/'.join(x))
        data = pd.concat([data['output_id'], data['index']]).sort_index(kind="mergesort").reset_index(drop = True)
        data = pd.DataFrame(data, columns=filepointer.background_peptide_fp['columns'])
        save_pd_toparquet(path, data, compression, verbose)


def save_forgrd_pep_dict(data, filepointer, compression=None, outbase=None, save_fasta=False, verbose=False):
    path_meta = switch_tmp_path(filepointer.junction_meta_fp, outbase)
    if data:
        data = pd.DataFrame(data.values(), index = data.keys()).reset_index()
        data.columns = filepointer.junction_meta_fp['columns']
        if save_fasta:
            path_fa = switch_tmp_path(filepointer.junction_peptide_fp, outbase)
            fasta = pd.concat([data['id'].apply(lambda x: '>' + '/'.join(x)), data['peptide']]).sort_index(kind="mergesort").reset_index(drop = True)
            fasta = pd.DataFrame(fasta, columns=filepointer.junction_peptide_fp['columns'])
            save_pd_toparquet(path_fa, fasta, compression, verbose)
            del fasta
        data = data.set_index('peptide')
        data = data.applymap(_convert_list_to_str)
        data = data.reset_index(drop = False)
        save_pd_toparquet(path_meta, data, compression, verbose)  # Test keep peptide name in the metadata file


def save_gene_expr_distr(data, filepointer, outbase, compression, verbose):
    if data:
        data = pd.DataFrame(data)
        path = switch_tmp_path(filepointer.gene_expr_fp, outbase)
        data.columns = filepointer.gene_expr_fp['columns']
        save_pd_toparquet(path, data, compression, verbose)


def save_kmer_matrix(data, graph_samples, filepointer, compression=None, outbase=None, verbose=False):
   '''The function saves the kmer segment and junction matrices per gene. The filepointer is kept open until all genes from the batch are written '''
   if data[0]:
       segm_path = switch_tmp_path(filepointer.kmer_segm_expr_fp, outbase)
       edge_path = switch_tmp_path(filepointer.kmer_edge_expr_fp, outbase)
       data[2] = pd.DataFrame(data[2], columns = graph_samples )
       data[3] = pd.DataFrame(data[3], columns = graph_samples)
       filepointer.kmer_segm_expr_fp['pqwriter'] = save_pd_toparquet(segm_path, pd.concat([pd.DataFrame({'kmer': data[0], 'is_cross_junction': data[1]}), data[2]], axis=1),
                         compression=compression, verbose=verbose, pqwriter=filepointer.kmer_segm_expr_fp['pqwriter'], writer_close=False)
       filepointer.kmer_edge_expr_fp['pqwriter'] = save_pd_toparquet(edge_path, pd.concat([pd.DataFrame({'kmer': data[0], 'is_cross_junction': data[1]}), data[3]], axis=1),
                         compression=compression, verbose=verbose, pqwriter=filepointer.kmer_edge_expr_fp['pqwriter'], writer_close=False)


def initialize_fp(output_path, mutation_mode, gzip_tag,
                  kmer_list, output_fasta, cross_graph_expr):

    ### Paths
    background_peptide_file_path = os.path.join(output_path, mutation_mode+ '_annot_peptides.fa.pq' + gzip_tag)
    background_kmer_file_path = os.path.join(output_path, mutation_mode + '_annot_kmer.pq' + gzip_tag)
    gene_expr_file_path = os.path.join(output_path, 'gene_expression_detail.pq' + gzip_tag)
    junction_meta_file_path = os.path.join(output_path, mutation_mode + '_sample_peptides_meta.tsv.gz.pq')
    junction_peptide_file_path = os.path.join(output_path, mutation_mode + '_sample_peptides.fa.pq' + gzip_tag)
    junction_kmer_file_path = os.path.join(output_path, mutation_mode + '_sample_kmer.pq' + gzip_tag)
    graph_kmer_segment_expr_path = os.path.join(output_path, mutation_mode + '_graph_kmer_SegmExpr.pq' + gzip_tag)
    graph_kmer_junction_expr_path = os.path.join(output_path, mutation_mode + '_graph_kmer_JuncExpr.pq' + gzip_tag)

    ### Fields
    fields_backgrd_pep_dict = ['fasta']
    fields_backgrd_kmer_dict = ['kmer']
    fields_gene_expr_file = ['gene', 'total_expr']
    fields_meta_peptide_dict = ['peptide', 'id', 'read_frame', 'gene_name', 'gene_chr', 'gene_strand',
                                'mutation_mode',
                                'junction_annotated', 'has_stop_codon', 'is_in_junction_list',
                                'is_isolated', 'variant_comb', 'variant_seg_expr', 'modified_exons_coord',
                                'original_exons_coord', 'vertex_idx',
                                'kmer_type']
    fields_forgrd_pep_dict = ['fasta']
    fields_forgrd_kmer_dict = ['kmer', 'id', 'segment_expr', 'is_cross_junction', 'junction_expr']
    fields_kmer_expr = ['kmer', 'is_cross_junction']
    fields_pept_expr = ['peptide', 'is_cross_junction']

    ### Grouping dict
    if output_fasta:    # Foreground peptide fasta - optional
        peptide_fp = output_info(junction_peptide_file_path, fields_forgrd_pep_dict)
    else:
        peptide_fp = None
    if cross_graph_expr:# Expression matrices with counts from full graph - optional
        junction_kmer_fp = None
        kmer_segm_expr_fp = output_info(graph_kmer_segment_expr_path, fields_kmer_expr, pq_writer=True)
        kmer_edge_expr_fp = output_info(graph_kmer_junction_expr_path, fields_kmer_expr, pq_writer=True)
    else: # Expression kmer information from single sample
        fields_meta_peptide_dict.insert(-1, 'junction_expr')
        fields_meta_peptide_dict.insert(-1, 'segment_expr')
        junction_kmer_fp = output_info(junction_kmer_file_path, fields_forgrd_kmer_dict, kmer_list )
        kmer_segm_expr_fp = None
        kmer_edge_expr_fp = None
    meta_peptide_fp = output_info(junction_meta_file_path, fields_meta_peptide_dict)
    background_fp = output_info(background_peptide_file_path, fields_backgrd_pep_dict )
    background_kmer_fp = output_info(background_kmer_file_path, fields_backgrd_kmer_dict, kmer_list )
    gene_expr_fp = output_info(gene_expr_file_path, fields_gene_expr_file)

    ### Filepointer object
    filepointer = Filepointer(peptide_fp, meta_peptide_fp, background_fp, junction_kmer_fp, background_kmer_fp,
                              gene_expr_fp, kmer_segm_expr_fp, kmer_edge_expr_fp)
    return filepointer


def output_info(path, file_columns, kmer_list=None, pq_writer=False):
    if kmer_list: # variable number of kmer files
        path_dict = {}
        for kmer_length in kmer_list:
            path_dict[kmer_length] = path.replace('kmer.pq', '{}mer.pq'.format(kmer_length))
        path = path_dict
    file_info = {'path': path, 'columns': file_columns}
    if pq_writer:
        file_info['pqwriter'] = None
    return file_info


def save_pd_toparquet(path, pd_df, compression = None, verbose = False, pqwriter=None, writer_close=True):
    s1 = timeit.default_timer()
    table = pa.Table.from_pandas(pd_df, preserve_index=False)

    if pqwriter is None:
        pqwriter = pq.ParquetWriter(path, table.schema, compression=compression)
    pqwriter.write_table(table)
    if writer_close:
        pqwriter.close()

    if verbose:
        file_name = os.path.basename(path)
        tot_shape = pd_df.shape[0]
        logging.info('Saving parquet {} with {} lines. Took {} seconds'.format(file_name, tot_shape,
                                                                      round(timeit.default_timer() - s1, 4)))
    return pqwriter

# def save_dict_toparquet(path, my_dict, columns, compression = None, verbose = False):
#     s1 = timeit.default_timer()
#     table = pa.Table.from_arrays([my_dict], columns["kmer"])
#     pa.parquet.write_table(table=my_dict, where=path, compression=compression, use_dictionary=columns)
#     if verbose:
#         file_name = os.path.basename(path)
#         tot_shape = len(my_dict[columns[0]])
#         logging.info('Saving parquet {} with {} lines. Took {} seconds'.format(file_name, tot_shape,
#                                                                       timeit.default_timer() - s1))


def collect_results(filepointer_item, outbase, compression, mutation_mode,  kmer_list = None):
    if filepointer_item is not None:
        s1 = timeit.default_timer()

        if kmer_list:
           file_to_collect = filepointer_item['path'].values()
        else:
            file_to_collect = [filepointer_item['path']]
        for file_path in file_to_collect:
            file_name = os.path.basename(file_path)
            tmp_file_list = glob.glob(os.path.join(outbase, 'tmp_out_{}_[0-9]*'.format(mutation_mode), file_name))
            tot_shape = 0
            for tmp_file in tmp_file_list:
                table = pq.read_table(tmp_file)
                if tot_shape == 0:
                    pqwriter = pq.ParquetWriter(file_path, table.schema, compression=compression)
                pqwriter.write_table(table)
                tot_shape += table.shape[0]
            if tmp_file_list:
                pqwriter.close()
                logging.info('Collecting {} with {} lines. Took {} seconds'.format(file_name, tot_shape, timeit.default_timer()-s1))


def remove_folder_list(commun_base):
    folder_list = glob.glob(commun_base + '*')
    for dir_path in folder_list:
        shutil.rmtree(dir_path, ignore_errors=True)


def read_pq_with_dict(file_path, columns):
    return pa.parquet.read_table(file_path, read_dictionary=columns)

def read_pq_with_pd(file_path):
    return pa.parquet.read_table(file_path).to_pandas()
