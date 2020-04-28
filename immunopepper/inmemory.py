# Function to be moved later to correct folders
from collections import OrderedDict
import os
import pandas as pd
import timeit

from .io_ import convert_namedtuple_to_str
from .io_ import _convert_list_to_str
from .utils import unpickler



def add_dict_kmer_forgrd(foregrd_dict, _namedtuple_list, filter_dict):

    for _namedtuple_kmer in _namedtuple_list:
        ord_dict = _namedtuple_kmer._asdict()
        del ord_dict['kmer']

        if _namedtuple_kmer.kmer in filter_dict:
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
        start_time = timeit.default_timer()
        backgrd_dict[_namedtuple_kmer.kmer] = 0
    return backgrd_dict


def filter_onkey_dict(dict_foregr, dict_back):
    pre_filt_kmers = dict_foregr.keys()
    for key_ in  pre_filt_kmers:
        if key_ in dict_back:
            del dict_foregr[key_]
    return dict_foregr

def add_dict_peptide(dict_peptides, _namedtuple_list ):
    for _namedtuple_peptide in _namedtuple_list:
        meta_data =  dict(_namedtuple_peptide._asdict())
        del meta_data['peptide']
        dict_peptides[_namedtuple_peptide.peptide] = meta_data
        return dict_peptides


def sort_records_kmer(nestedlist_of_namedtuples):
    records = dict([(record.kmer, record) for transcript_list in nestedlist_of_namedtuples for record in
                    transcript_list])
    records = OrderedDict(sorted(records.items()))
    return records



def write_gene_result(gene_result, dict_pept_forgrd, dict_pept_backgrd, dict_kmer_foregr, dict_kmer_back, logging):

    if len(gene_result['background_peptide_list']):
        records = gene_result['background_peptide_list']
        dict_pept_backgrd = add_dict_peptide(dict_pept_backgrd, records)

    if len(gene_result['background_kmer_lists']):
        records = [item for sublist in gene_result['background_kmer_lists'] for item in sublist]
        dict_kmer_back = add_dict_kmer_back(dict_kmer_back, records)

    if len(gene_result['output_metadata_list']):
        records = gene_result['output_metadata_list']
        dict_pept_forgrd = add_dict_peptide(dict_pept_forgrd, records)

    if len(gene_result['output_kmer_lists']):
        records = [item for sublist in gene_result['output_kmer_lists'] for item in sublist]
        dict_kmer_foregr = add_dict_kmer_forgrd(dict_kmer_foregr, records,  dict_kmer_back)


    return dict_pept_forgrd, dict_pept_backgrd, dict_kmer_foregr, dict_kmer_back



def save_backgrd_kmer_dict(dict, save_path, compression = None):
        df = pd.DataFrame(dict.keys(), columns = ['kmer'])
        df.to_parquet(save_path,  engine='fastparquet',
                  compression=compression, index=False)

def save_forgrd_kmer_dict(dict, save_path, compression = None):
        df = pd.DataFrame(dict.values(), index = dict.keys())
        df = df.applymap(repr)
        df = df.rename_axis('kmer').reset_index()
        df.to_parquet(save_path, engine='fastparquet',
                  compression=compression, index=False)

def save_backgrd_pep_dict(dict, save_path_back_pep, compression = None):
    fasta = pd.DataFrame(dict.values(), index = dict.keys()).reset_index()
    fasta = pd.concat([fasta['id'], fasta['index']]).sort_index().reset_index(drop = True)
    fasta = pd.DataFrame(fasta, columns=['fasta'])
    fasta.to_parquet(save_path_back_pep, engine='fastparquet',
                  compression=compression, index=False)

def save_forgrd_pep_dict(dict, save_path_forgr_pep, save_path_meta_pep, compression = None):
    df = pd.DataFrame(dict.values(), index = dict.keys()).reset_index()
    fasta = pd.concat([df['output_id'], df['index']]).sort_index().reset_index(drop = True)
    fasta = pd.DataFrame(fasta, columns=['fasta'])
    df = df.drop(['output_id', 'index'], axis=1)
    fasta.to_parquet(save_path_forgr_pep, engine='fastparquet',
                  compression=compression, index=False)
    del fasta
    df['original_exons_coord'] = df['original_exons_coord'].apply(convert_namedtuple_to_str, args=(None, ';'))
    df['modified_exons_coord'] = df['modified_exons_coord'].apply(convert_namedtuple_to_str, args=(None, ';'))
    df = df.applymap(repr)
    df.to_parquet(save_path_meta_pep, engine='fastparquet',
                  compression=compression, index=False)
