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


def save_backgrd_kmer_set(set_, filepointer, kmer_length, compression=None, outbase=None, verbose=False):
    path = switch_tmp_path(filepointer.background_kmer_fp, outbase, kmer_length)
    if set_:
        df = pd.DataFrame(set_, columns=['kmer'])
        save_pd_toparquet(path, df, compression, verbose)

def save_forgrd_kmer_dict(dict_, filepointer, kmer_length, compression=None, outbase=None, verbose=False):
    path = switch_tmp_path(filepointer.junction_kmer_fp, outbase, kmer_length)
    if dict_:
        df = pd.DataFrame(dict_.values(), index=dict_.keys())
        df = df.applymap(_convert_list_to_str)
        df = df.rename_axis('kmer').reset_index()
        df = df.loc[:, filepointer.junction_kmer_fp['columns']]
        save_pd_toparquet(path, df, compression, verbose)

def save_backgrd_pep_dict(dict_, filepointer, compression = None, outbase=None, verbose=False):
    path = switch_tmp_path(filepointer.background_peptide_fp, outbase)
    if dict_:
        fasta = pd.DataFrame(dict_.values(), index = dict_.keys()).reset_index()
        fasta['output_id'] = fasta['output_id'].apply(lambda x: '>' + '/'.join(x))
        fasta = pd.concat([fasta['output_id'], fasta['index']]).sort_index(kind="mergesort").reset_index(drop = True)
        fasta = pd.DataFrame(fasta, columns=['fasta'])
        save_pd_toparquet(path, fasta, compression, verbose)

def save_forgrd_pep_dict(dict_, filepointer, compression=None, outbase=None, verbose=False):
    path_fa = switch_tmp_path(filepointer.junction_peptide_fp, outbase)
    path_meta = switch_tmp_path(filepointer.junction_meta_fp, outbase)
    if dict_:
        df = pd.DataFrame(dict_.values(), index = dict_.keys()).reset_index()
        df['output_id'] = df['output_id'].apply(lambda x: '>' + '/'.join(x))
        fasta = pd.concat([df['output_id'], df['index']]).sort_index(kind="mergesort").reset_index(drop = True)
        fasta = pd.DataFrame(fasta, columns=['fasta'])
        save_pd_toparquet(path_fa, fasta, compression, verbose)
        del fasta
        df = df.drop(['output_id', 'index'], axis=1)
        df = df.loc[:, filepointer.junction_meta_fp['columns']]
        df = df.applymap(_convert_list_to_str)
        save_pd_toparquet(path_meta, df, compression, verbose)

def save_gene_expr_distr(gene_expr_distr_list, filepointer, outbase, compression, verbose):
    df = pd.DataFrame(gene_expr_distr_list)
    path = switch_tmp_path(filepointer.gene_expr_fp, outbase)
    df.columns = filepointer.gene_expr_fp['columns']
    save_pd_toparquet(path, df, compression, verbose)


def initialize_fp(junction_peptide_file_path, junction_meta_file_path, background_peptide_file_path,
                    junction_kmer_file_path, background_kmer_file_path, gene_expr_file_path, kmer_list):

    fields_forgrd_pep_dict = ['fasta']
    fields_meta_peptide_dict = ['id','read_frame','gene_name','gene_chr','gene_strand','mutation_mode',
                                'peptide_annotated','junction_annotated','has_stop_codon','is_in_junction_list',
                                'is_isolated','variant_comb','variant_seg_expr','modified_exons_coord',
                                'original_exons_coord','vertex_idx','junction_expr','segment_expr', 'kmer_type']
    fields_backgrd_pep_dict = ['fasta']
    fields_forgrd_kmer_dict = ['kmer', 'id', 'expr', 'is_cross_junction', 'junction_count']
    fields_backgrd_kmer_dict = ['kmer']
    fields_gene_expr_file = ['gene', 'total_expr']


    peptide_fp = output_info(junction_peptide_file_path, fields_forgrd_pep_dict)
    meta_peptide_fp = output_info(junction_meta_file_path, fields_meta_peptide_dict )
    background_fp = output_info(background_peptide_file_path, fields_backgrd_pep_dict )
    junction_kmer_fp = output_info(junction_kmer_file_path,fields_forgrd_kmer_dict , kmer_list )
    background_kmer_fp = output_info(background_kmer_file_path, fields_backgrd_kmer_dict ,kmer_list )
    gene_expr_fp = output_info(gene_expr_file_path, fields_gene_expr_file)
    filepointer = Filepointer(peptide_fp, meta_peptide_fp, background_fp, junction_kmer_fp, background_kmer_fp, gene_expr_fp )
    return filepointer


def output_info(path, file_columns, kmer_list=None):
    if kmer_list: # variable number of kmer files
        path_dict = {}
        for kmer_length in kmer_list:
            path_dict[kmer_length] = path.replace('kmer.pq', '{}mer.pq'.format(kmer_length))
        path = path_dict
    file_info = {'path': path, 'columns': file_columns}
    return file_info


def save_pd_toparquet(path, pd_df, compression = None, verbose = False):
    s1 = timeit.default_timer()
    table = pa.Table.from_pandas(pd_df, preserve_index=False)
    pqwriter = pq.ParquetWriter(path, table.schema, compression=compression)
    pqwriter.write_table(table)
    pqwriter.close()
    if verbose:
        file_name = os.path.basename(path)
        tot_shape = pd_df.shape[0]
        logging.info('Saving parquet {} with {} lines. Took {} seconds'.format(file_name, tot_shape,
                                                                      timeit.default_timer() - s1))


def collect_results(filepointer_item, outbase, compression, mutation_mode,  kmer_list = None):
    s1 = timeit.default_timer()

    if kmer_list:
       file_to_collect = filepointer_item['path'].values()
    else:
        file_to_collect = [filepointer_item['path']]


    for file_path in file_to_collect:
        file_name = os.path.basename(file_path)
        tmp_file_list = glob.glob(os.path.join(outbase, 'tmp_out_{}*'.format(mutation_mode), file_name))
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
