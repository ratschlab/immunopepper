# Function to be moved later to correct folders

import os
import pandas as pd

import datrie
from .io import convert_namedtuple_to_str
from .io import _convert_list_to_str
from .utils import unpickler

def create_kmer_trie(base = False):
    if base:
        trie = datrie.BaseTrie(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
    else:
        trie = datrie.Trie(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
    return trie


def add_trie_kmer_forgrd(trie, _namedtuple_list, kmer_field_list, filter_trie):
    base_dict = {}
    for _namedtuple_kmer in _namedtuple_list:
        meta_data = convert_namedtuple_to_str(_namedtuple_kmer, kmer_field_list[1:])
        if _namedtuple_kmer.kmer in filter_trie: #TODO add Back
            continue
        if _namedtuple_kmer.kmer not in trie: #potential slow down here
            for field in kmer_field_list[1:]:
                base_dict[field] = [ getattr(_namedtuple_kmer, field) ]
            trie[_namedtuple_kmer.kmer] = base_dict
        else:
            for field in kmer_field_list[1:]:
                item_= getattr(_namedtuple_kmer, field)
                if item_ not in trie[_namedtuple_kmer.kmer][field]:
                    trie[_namedtuple_kmer.kmer][field].append(item_)
    return trie

def add_trie_kmer_back(trie, _namedtuple_list):
    for _namedtuple_kmer in _namedtuple_list:
        trie[_namedtuple_kmer.kmer] = 0
    return trie

def filter_onkey_trie(trie_foregr, trie_back):
    for key_ in  trie_foregr:
        if key_ in trie_back:
            del trie_foregr[key_]
    return trie_foregr

def add_trie_peptide(trie, _namedtuple ):
    meta_data =  dict(_namedtuple._asdict())
    del meta_data['peptide']
    trie[_namedtuple.peptide] = meta_data
    return trie


def write_gene_result(gene_result, trie_pept_forgrd, trie_pept_backgrd, trie_kmer_foregr, trie_kmer_back):

    ### define fields relevant for output
    #back_pep_field_list = ['id', 'new_line', 'peptide']
    for peptide in gene_result['background_peptide_list']:
        trie_pept_backgrd = add_trie_peptide(trie_pept_backgrd, peptide)

    #back_kmer_field_list = ['kmer', 'id', 'expr', 'is_cross_junction']
    for kmer_list in gene_result['background_kmer_lists']:
        trie_kmer_back = add_trie_kmer_back(trie_kmer_back, kmer_list)

    #junc_pep_field_list = ['output_id', 'id', 'new_line', 'peptide']
<<<<<<< HEAD
    if os.path.exists(gene_result['output_metadata_list']):
        with open(gene_result['output_metadata_list'], 'rb') as fh:
            for record in unpickler(fh):
                trie_pept_forgrd = add_trie_peptide(trie_pept_forgrd, record)
        os.remove(gene_result['output_metadata_list'])
=======
        for record in gene_result['output_metadata_list']:
            trie_pept_forgrd = add_trie_peptide(trie_pept_forgrd, record)
>>>>>>> 3a408cb82f611351784fd20289f3e4d7b486b240


        kmer_field_list = ['kmer', 'id', 'expr', 'is_cross_junction', 'junction_count']
        for kmer_list in gene_result['output_kmer_lists']:
            trie_kmer_foregr = add_trie_kmer_forgrd(trie_kmer_foregr, kmer_list, kmer_field_list,
                                                            trie_kmer_back)


    return trie_pept_forgrd, trie_pept_backgrd, trie_kmer_foregr, trie_kmer_back



def save_backgrd_kmer_trie(trie, save_path, compression = None):
        df = pd.DataFrame(trie.keys(), columns = ['kmer'])
        df.to_parquet(save_path,  engine='fastparquet',
                  compression=compression, index=False)

def save_forgrd_kmer_trie(trie, save_path, compression = None):
        df = pd.DataFrame(trie.values(), index = trie.keys())
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
                  compression=compression)
    del fasta
    # TODO potential slow down here?
    df['original_exons_coord'] = df['original_exons_coord'].apply(convert_namedtuple_to_str, args=(None, ';'))
    df['modified_exons_coord'] = df['modified_exons_coord'].apply(convert_namedtuple_to_str, args=(None, ';'))
    df['junction_expr'] = df['junction_expr'].apply(_convert_list_to_str)
    df['variant_comb'] = df['variant_comb'].apply(_convert_list_to_str)
    df['variant_seg_expr'] = df['variant_seg_expr'].apply(_convert_list_to_str)
    df['vertex_idx'] = df['vertex_idx'].apply(_convert_list_to_str)
    df['junction_expr'] = df['junction_expr'].apply(_convert_list_to_str)
    df.to_parquet(save_path_meta_pep, engine='fastparquet',
                  compression=compression)
    del df