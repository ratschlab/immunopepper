# Function to be moved later to correct folders
import glob
import os
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import shutil
import timeit

from .io_ import convert_namedtuple_to_str
from .namedtuples import Filepointer

### Move to Filter
def add_dict_kmer_forgrd(foregrd_dict, _namedtuple_list, filter_dict, remove_annot=True):

    for _namedtuple_kmer in _namedtuple_list:
        ### Prepare metadata
        ord_dict = _namedtuple_kmer._asdict()
        del ord_dict['kmer']

        ### Remove annotaion
        if (remove_annot) and (_namedtuple_kmer.kmer) in filter_dict:
            continue

        if _namedtuple_kmer.kmer not in foregrd_dict:
            dic_with_sets = dict(zip(ord_dict.keys(), [{i} for i in ord_dict.values()]))
            foregrd_dict[_namedtuple_kmer.kmer] = dic_with_sets
        else:
            for field, value in ord_dict.items():
                foregrd_dict[_namedtuple_kmer.kmer][field].add(value)

    return foregrd_dict


def add_dict_kmer_back(backgrd_dict, _namedtuple_list):
    for _namedtuple_kmer in _namedtuple_list:
        if _namedtuple_kmer.kmer in backgrd_dict:
            continue
        backgrd_dict[_namedtuple_kmer.kmer] = 0
    return backgrd_dict


def add_dict_peptide(dict_peptides, _namedtuple_list ):
    for _namedtuple_peptide in _namedtuple_list:
        meta_data =  dict(_namedtuple_peptide._asdict())
        del meta_data['peptide']
        dict_peptides[_namedtuple_peptide.peptide] = meta_data
        return dict_peptides


def filter_onkey_dict(dict_foregr, dict_back):
    pre_filt_kmers = list(dict_foregr.keys())
    for key_ in  pre_filt_kmers:
        if key_ in dict_back:
            del dict_foregr[key_]
    return dict_foregr

### Move to mode_build
def write_gene_result(gene_result, dict_pept_forgrd, dict_pept_backgrd, dict_kmer_foregr, dict_kmer_back, filepointer,
                      remove_annot = True, outbase = None):


    if len(gene_result['background_peptide_list']):
        records = gene_result['background_peptide_list']
        dict_pept_backgrd = add_dict_peptide(dict_pept_backgrd, records)
        save_backgrd_pep_dict(dict_pept_backgrd, filepointer, compression=None, outbase = outbase)
        dict_pept_backgrd = {}

    if len(gene_result['background_kmer_lists']):
        records = [item for sublist in gene_result['background_kmer_lists'] for item in sublist]
        dict_kmer_back = add_dict_kmer_back(dict_kmer_back, records)
        if not remove_annot:
            save_backgrd_kmer_dict(dict_kmer_back, filepointer, compression=None, outbase = outbase)
            dict_kmer_back = {}


    if len(gene_result['output_metadata_list']):
        records = gene_result['output_metadata_list']
        dict_pept_forgrd = add_dict_peptide(dict_pept_forgrd, records)
        if not remove_annot:
            save_forgrd_pep_dict(dict_pept_forgrd, filepointer, compression=None, outbase = outbase)
            dict_pept_forgrd = {}

    if len(gene_result['output_kmer_lists']):
        records = [item for sublist in gene_result['output_kmer_lists'] for item in sublist]
        dict_kmer_foregr = add_dict_kmer_forgrd(dict_kmer_foregr, records,  dict_kmer_back, remove_annot)
        if not remove_annot:
            save_forgrd_kmer_dict(dict_kmer_foregr, filepointer, compression=None, outbase = outbase)
            dict_kmer_foregr = {}

    return dict_pept_forgrd, dict_pept_backgrd, dict_kmer_foregr, dict_kmer_back

### Move to io_
def switch_tmp_path(filepointer_item, outbase=None):
    if outbase:
        fp = None
        path = os.path.join(outbase, os.path.basename(filepointer_item['path']))
    else:
        fp = filepointer_item['filepointer']
        path = filepointer_item['path']
    return fp, path


def save_backgrd_kmer_dict(dict_, filepointer, compression = None, outbase = None):
    fp, path = switch_tmp_path(filepointer.background_kmer_fp, outbase)
    if dict_:
        df = pd.DataFrame(dict_.keys(), columns = ['kmer'])
        save_pd_toparquet(fp, path, df, compression)

def save_forgrd_kmer_dict(dict_, filepointer, compression = None, outbase=None):
    fp, path = switch_tmp_path(filepointer.junction_kmer_fp, outbase)
    if dict_:
        df = pd.DataFrame(dict_.values(), index = dict_.keys())
        df = df.applymap(repr)
        df = df.rename_axis('kmer').reset_index()
        if filepointer.junction_kmer_fp['columns']:
            df = df.loc[:, filepointer.junction_kmer_fp['columns'] ]
        save_pd_toparquet(fp, path, df, compression)

def save_backgrd_pep_dict(dict_, filepointer, compression = None, outbase=None):
    fp, path = switch_tmp_path(filepointer.background_peptide_fp, outbase)
    if dict_:
        fasta = pd.DataFrame(dict_.values(), index = dict_.keys()).reset_index()
        fasta = pd.concat([fasta['output_id'], fasta['index']]).sort_index().reset_index(drop = True)
        fasta = pd.DataFrame(fasta, columns=['fasta'])
        save_pd_toparquet(fp, path, fasta, compression)

def save_forgrd_pep_dict(dict_, filepointer, compression = None, outbase=None):
    fp_fa, path_fa = switch_tmp_path(filepointer.junction_peptide_fp, outbase)
    fp_meta, path_meta = switch_tmp_path(filepointer.junction_meta_fp, outbase)
    if dict_:
        df = pd.DataFrame(dict_.values(), index = dict_.keys()).reset_index()
        fasta = pd.concat([df['output_id'], df['index']]).sort_index().reset_index(drop = True)
        fasta = pd.DataFrame(fasta, columns=['fasta'])
        df = df.drop(['output_id', 'index'], axis=1)
        save_pd_toparquet(fp_fa, path_fa, fasta, compression)
        del fasta
        df['original_exons_coord'] = df['original_exons_coord'].apply(convert_namedtuple_to_str, args=(None, ';'))
        df['modified_exons_coord'] = df['modified_exons_coord'].apply(convert_namedtuple_to_str, args=(None, ';'))
        df = df.applymap(repr)
        if filepointer.junction_meta_fp['columns']:
            df = df.loc[:, filepointer.junction_meta_fp['columns'] ]
        save_pd_toparquet(fp_meta, path_meta, df, compression)


def initialize_parquet(junction_peptide_file_path, junction_meta_file_path, background_peptide_file_path,
                    junction_kmer_file_path, background_kmer_file_path , remove_annot=True):
    fields_forgrd_pep_dict = ['fasta']
    fields_meta_peptide_dict = ['id','read_frame','gene_name','gene_chr','gene_strand','mutation_mode',
                                'peptide_annotated','junction_annotated','has_stop_codon','is_in_junction_list',
                                'is_isolated','variant_comb','variant_seg_expr','modified_exons_coord',
                                'original_exons_coord','vertex_idx','junction_expr','segment_expr']
    fields_backgrd_pep_dict = ['fasta']
    fields_forgrd_kmer_dict = ['kmer', 'id', 'expr', 'is_cross_junction', 'junction_count'] ## Needs to be in memory, Unless no filtering
    fields_backgrd_kmer_dict = ['kmer']

    peptide_fp = fp_with_pq_schema(junction_peptide_file_path, None)#fields_forgrd_pep_dict ) #TODO add back file pointer creation
    meta_peptide_fp = fp_with_pq_schema(junction_meta_file_path, None) # fields_meta_peptide_dict)
    background_fp = fp_with_pq_schema(background_peptide_file_path, None) #fields_backgrd_pep_dict)
    junction_kmer_fp = fp_with_pq_schema(junction_kmer_file_path, None ) #fields_forgrd_kmer_dict)
    background_kmer_fp = fp_with_pq_schema(background_kmer_file_path, None) #fields_backgrd_kmer_dict)
    filepointer = Filepointer(peptide_fp, meta_peptide_fp, background_fp, junction_kmer_fp, background_kmer_fp)
    return filepointer


def fp_with_pq_schema(path, file_columns, compression = None):
    if file_columns:
        data = {col: repr(set({})) for col in file_columns}
        data = pd.DataFrame(data, index=[0])
        table = pa.Table.from_pandas(data,  preserve_index=False)
        pqwriter = pq.ParquetWriter(path, table.schema, compression)
        file_info = {'path':path,'filepointer':pqwriter, 'columns':file_columns}
    else:
        file_info = {'path': path,'filepointer': None, 'columns': None}
        print('file pointer is not initialized')
    return file_info


def save_pd_toparquet(filepointer, path, pd_df, compression = None):
    table = pa.Table.from_pandas(pd_df, preserve_index=False)
    if filepointer is None:
        pqwriter = pq.ParquetWriter(path, table.schema, compression)
        pqwriter.write_table(table)
        pqwriter.close()
    else:
        pqwriter = filepointer
        pqwriter.write_table(table)

def collect_results(filepointer_item,outbase,logging):
    s1 = timeit.default_timer()
    file_name = os.path.basename(filepointer_item['path'])
    tmp_file_list = glob.glob(os.path.join(outbase,'tmp_out_*',file_name))
    pqwriter = filepointer_item['filepointer']
    tot_shape = 0
    for tmp_file in tmp_file_list:
        table = pq.read_table(tmp_file)
        if tot_shape == 0: #TODO do not keep open filepointers
            pqwriter = pq.ParquetWriter(file_name, table.schema, compression = None)
        pqwriter.write_table(table)
        tot_shape += table.shape[0]
    if tmp_file_list:
        logging.info('Collecting {} with {} lines. Took {} seconds'.format(file_name, tot_shape, timeit.default_timer()-s1))


def remove_folder_list(commun_base):
    folder_list = glob.glob(commun_base + '*')
    for dir_path in folder_list:
        shutil.rmtree(dir_path, ignore_errors=True)
