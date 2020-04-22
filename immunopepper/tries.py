# Function to be moved later to correct folders
from collections import OrderedDict
import datrie
import os
import pandas as pd
import timeit

from .io import convert_namedtuple_to_str
from .io import _convert_list_to_str
from .utils import unpickler


def create_kmer_trie(base = False):
    if base:
        trie = datrie.BaseTrie(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
    else:
        trie = datrie.Trie(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
    return trie


def add_trie_kmer_forgrd(trie, _namedtuple_od, filter_trie):

    for _namedtuple_kmer in _namedtuple_od.values():
        ord_dict = _namedtuple_kmer._asdict()
        del ord_dict['kmer']

        if _namedtuple_kmer.kmer in filter_trie: #TODO add Back
            continue
        if _namedtuple_kmer.kmer not in trie: #potential slow down here

            dic_with_sets = dict(zip(ord_dict.keys(), [{i} for i in ord_dict.values()]))
            trie[_namedtuple_kmer.kmer] = dic_with_sets
        else:
            for field, value in ord_dict.items():
                trie[_namedtuple_kmer.kmer][field].add(value)

    return trie


def add_trie_kmer_back(trie, _namedtuple_od, logging):
    for _namedtuple_kmer in _namedtuple_od.values():
        if _namedtuple_kmer.kmer in trie: 
            continue
        start_time = timeit.default_timer()
        trie[_namedtuple_kmer.kmer] = 0
        logging.info('background_kmer trie single  update {}'.format(timeit.default_timer() - start_time))
    return trie


def filter_onkey_trie(trie_foregr, trie_back):
    for key_ in  trie_foregr:
        if key_ in trie_back:
            del trie_foregr[key_]
    return trie_foregr

def add_trie_peptide(trie, _namedtuple_od ):
    for _namedtuple_peptide in _namedtuple_od.values():
        meta_data =  dict(_namedtuple_peptide._asdict())
        del meta_data['peptide']
        trie[_namedtuple_peptide.peptide] = meta_data
        return trie


def sort_records_kmer(nestedlist_of_namedtuples):
    records = dict([(record.kmer, record) for transcript_list in nestedlist_of_namedtuples for record in
                    transcript_list])
    records = OrderedDict(sorted(records.items()))
    return records

def sort_records_peptide(list_of_namedtuples):
    records = dict([(record.peptide, record) for record in list_of_namedtuples])
    records = OrderedDict(sorted(records.items()))
    return records


def write_gene_result(gene_result, trie_pept_forgrd, trie_pept_backgrd, trie_kmer_foregr, trie_kmer_back, logging):

    if len(gene_result['background_peptide_list']):
        sorted_record_odict = sort_records_peptide(gene_result['background_peptide_list'])
        trie_pept_backgrd = add_trie_peptide(trie_pept_backgrd, sorted_record_odict)

    if len(gene_result['background_kmer_lists']):
        sorted_record_odict = sort_records_kmer(gene_result['background_kmer_lists'])
        trie_kmer_back = add_trie_kmer_back(trie_kmer_back, sorted_record_odict, logging)

    if len(gene_result['output_metadata_list']):
        sorted_record_odict = sort_records_peptide(gene_result['output_metadata_list'])
        trie_pept_forgrd = add_trie_peptide(trie_pept_forgrd, sorted_record_odict)

    if len(gene_result['output_kmer_lists']):
        sorted_record_odict = sort_records_kmer(gene_result['output_kmer_lists'])
        trie_kmer_foregr = add_trie_kmer_forgrd(trie_kmer_foregr, sorted_record_odict,  trie_kmer_back)


    return trie_pept_forgrd, trie_pept_backgrd, trie_kmer_foregr, trie_kmer_back



def save_backgrd_kmer_trie(trie, save_path, compression = None):
        df = pd.DataFrame(trie.keys(), columns = ['kmer'])
        df.to_parquet(save_path,  engine='fastparquet',
                  compression=compression, index=False)

def save_forgrd_kmer_trie(trie, save_path, compression = None):
        df = pd.DataFrame(trie.values(), index = trie.keys())
        df = df.applymap(repr)
        df = df.rename_axis('kmer').reset_index()
        df.to_parquet(save_path, engine='fastparquet',
                  compression=compression, index=False)

def save_backgrd_pep_trie(trie, save_path_back_pep, compression = None):
    fasta = pd.DataFrame(trie.values(), index = trie.keys()).reset_index()
    fasta = pd.concat([fasta['id'], fasta['index']]).sort_index().reset_index(drop = True)
    fasta = pd.DataFrame(fasta, columns=['fasta'])
    fasta.to_parquet(save_path_back_pep, engine='fastparquet',
                  compression=compression, index=False)

def save_forgrd_pep_trie(trie, save_path_forgr_pep, save_path_meta_pep, compression = None):
    df = pd.DataFrame(trie.values(), index = trie.keys()).reset_index()
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
