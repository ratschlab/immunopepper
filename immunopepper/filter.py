"""Contains the functionality used for filtering"""
import sys
import re

import numpy as np

from .constant import NOT_EXIST

def _collect_remove_ids(exon_junction_dict):
    """ For all exon pairs around the same intron, keep the longest one.

    Parameters
    ----------
    exon_junction_dict: dict. (read_frame, pos_mid_1, pos_mid_2, variant_combination)
                               -> List[Tuple(output_idx, pos_start, pos_end,)]

    Returns
    -------
    remove_id_list: List[junction ids to remove]
    """

    # TODO: [Warning!] if two output lines have identical
    # coordinates and readframe both of them will be removed
    remove_id_list = []
    removed = set()
    for exon_pair_list in exon_junction_dict.values():
        L = len(exon_pair_list)
        if L < 2:
            continue
        for i, exon_pair in enumerate(exon_pair_list):
            if exon_pair[0] in removed:
                continue
            for j in range(i + 1, L):
                if exon_pair_list[j][0] in removed:
                    continue
                i_pos1 = exon_pair[1]
                i_pos2 = exon_pair[2]
                j_pos1 = exon_pair_list[j][1]
                j_pos2 = exon_pair_list[j][2]
                if i_pos1 >= j_pos1 and i_pos2 <= j_pos2:
                    remove_id_list.append(exon_pair[0])
                    removed.add(exon_pair[0])
                    break
                if j_pos1 >= i_pos1 and j_pos2 <= i_pos2:
                    remove_id_list.append(exon_pair_list[j][0])
                    removed.add(exon_pair_list[j][0])

    return remove_id_list


def _get_exon_junction_dict(metadata_list, strand):
    """ Creates an auxiliary dictionary used for filtering. Assigns to each junction all overlapping exon pairs.

    Parameters
    ----------
    metadata_list: List(Output_metadata).

    Returns
    -------
    exon_dict: dict. (read_frame, pos_mid_1, pos_mid_2, variant_combination)
        -> List[Tuple(output_idx, pos_start, pos_end,)]
    """
    exon_dict = {}
    for metadata in metadata_list:
        ## TODO: need to come up with a new way to index the exon dict
        coord = metadata.modified_exons_coord
        idx = metadata.output_id
        read_frame = metadata.read_frame.read_phase
        ### figure out intron vs other coordinates
        if strand == '+':
            i1, i2, o1, o2 = coord.stop_v1, coord.start_v2, coord.start_v1, coord.stop_v2
        else:
            i1, i2, o1, o2 = coord.stop_v2, coord.start_v1, coord.start_v2, coord.stop_v1
        key = (read_frame, i1, i2)
        if key in exon_dict:
            exon_dict[key].append((idx, o1, o2))
        else:
            exon_dict[key] = [(idx, o1, o2)]

    return exon_dict


def get_filtered_metadata_list(metadata_list, strand):
    """ Given a lis of exon junctions, remove the ones that redundantly cover a junction

    Parameters
    ----------
    metadata_list: List(Output_metadata),
    strand: strand of the gene

    Returns
    -------
    filtered_meetadata_list: List of metadata objects remaining after filter
    """

    exon_dict = _get_exon_junction_dict(metadata_list, strand)
    remove_id_list = _collect_remove_ids(exon_dict)
    return list(filter(lambda m: m.output_id not in remove_id_list, metadata_list))


def junction_is_annotated(gene, gene_to_transcript_table, transcript_to_cds_table):
    """ Indicates whether a junction also appears in any transcript given by .gtf file

    Parameters
    ----------
    gene: Object. Created by SplAdder.
    gene_to_transcript_table: Dict. "Gene" -> "transcript"
    transcript_to_cds_table: Dict. "transcript" -> "cds"

    Returns
    -------
    junction_flag: 2D array matrix with shape len(edges)*len(edges).
        True indicates the exon pair also appear in gtf file
        while False indicates the exon pair not exist or it does not show in the
        given .gtf file.
    """
    ts_list = gene_to_transcript_table[gene.name]
    vertices = gene.splicegraph.vertices
    edges = gene.splicegraph.edges

    junction_flag = np.zeros(edges.shape, dtype='bool')
    for ts in ts_list:
        if not ts in transcript_to_cds_table:
            continue

        curr_ts = transcript_to_cds_table[ts]
        if len(curr_ts) < 2:
            continue
        curr_ts = np.array(curr_ts)[:, [0, 1]]
        transcript_junctions = [':'.join(x) for x in
                                curr_ts.ravel()[1:-1].reshape(curr_ts.shape[0] - 1, 2).astype('str')]

        for x in np.array(np.where(np.triu(edges))).T:
            if '%i:%i' % (vertices[1, x[0]], vertices[0, x[1]]) in transcript_junctions:
                junction_flag[x[0], x[1]] = 1
                junction_flag[x[1], x[0]] = 1

    return junction_flag


def junction_tuple_is_annotated(junction_flag, vertex_id_tuple):
    """Returns 1 if any junction in the given tuple of vertexes is flagged and 0 otherwise.

    Parameters
    ----------
    junction_flag: 2-d array with shape (number_of_vertices, number_of_vertices). if (i,j) entry is 1, i-th vertex and
        j-th vertex are connected.
    vertex_id_tuple: triple-elements tuple. (v1,v2,v3) represents the vertex id for the junction pair.

    Returns
    -------
    int or list of int depends on it being a 2-vertex junction of 3-vertex junction.

    """
    if NOT_EXIST in vertex_id_tuple:
        return np.nan
    for i in range(len(vertex_id_tuple) - 1):
        if junction_flag[vertex_id_tuple[i], vertex_id_tuple[i+1]]:
            return 1
    return 0


def junction_is_in_given_list(splicegraph, vertex_id_list, strand, junction_list):
    """Check if the intron is in the user provided list of junctions"""

    if NOT_EXIST in vertex_id_list or junction_list is None:
        return NOT_EXIST

    vertex_coord_list = splicegraph.vertices
    if strand == '-':
        vertex_id_list = vertex_id_list[::-1] # in negative strand case, the vertex id list is [5,3] instead of [3,5]
    for i in range(len(vertex_id_list) - 1):
        junc_comb = ':'.join([str(vertex_coord_list[1, vertex_id_list[i]]),  str(vertex_coord_list[0, vertex_id_list[i+1]]), strand])
        if junc_comb in junction_list:
            return 1
    return 0


def peptide_is_annotated(background_peptide_list, peptide):
    """Check if the translated exon-pair peptide also appears in the background peptide translated from annotation file.

    Parameters
    ----------
    background_peptide_list: List(OutputBackground) All the peptide translated from transcripts in annotation file
    peptide: str. peptide translated from certain exon-pairs

    Returns
    -------
    count: int. Count how many matches exist between background peptide and given peptide
    """
    match_ts_list = []
    for background_peptide in background_peptide_list:
        ref_peptide = background_peptide.peptide
        transcript_id = background_peptide.output_id
        if not re.search(peptide, ref_peptide) is None:
            match_ts_list.append(transcript_id)
    return match_ts_list




def add_dict_kmer_forgrd(foregrd_dict, _namedtuple_list, filter_dict, remove_annot=True):
    """...

        Parameters
        ----------
        foregrd_dict:
        _namedtuple_list:
        filter_dict:

        """
    for _namedtuple_kmer in _namedtuple_list:
        ### Prepare metadata
        ord_dict = _namedtuple_kmer._asdict()
        del ord_dict['kmer']

        ### Remove annotation
        if (remove_annot) and (_namedtuple_kmer.kmer) in filter_dict:
            continue

        if _namedtuple_kmer.kmer not in foregrd_dict:
            dic_with_sets = dict(zip(ord_dict.keys(), [{i} for i in ord_dict.values()]))
            foregrd_dict[_namedtuple_kmer.kmer] = dic_with_sets
        else:
            for field, value in ord_dict.items():
                foregrd_dict[_namedtuple_kmer.kmer][field].add(value)



def add_dict_kmer_back(backgrd_dict, _namedtuple_list):
    """...

        Parameters
        ----------
        backgrd_dict:
        _namedtuple_list:

        """
    for _namedtuple_kmer in _namedtuple_list:
        if _namedtuple_kmer.kmer in backgrd_dict:
            continue
        backgrd_dict[_namedtuple_kmer.kmer] = 0


def add_dict_peptide(dict_peptides, _namedtuple_list ):
    """...

        Parameters
        ----------
        dict_peptides:
        _namedtuple_list:

        """
    for _namedtuple_peptide in _namedtuple_list:
        meta_data =  dict(_namedtuple_peptide._asdict())
        del meta_data['peptide']
        dict_peptides[_namedtuple_peptide.peptide] = meta_data


def filter_onkey_dict(dict_foregr, dict_back):
    """...

        Parameters
        ----------
        dict_foregr:
        dict_back:

        """
    pre_filt_kmers = list(dict_foregr.keys())
    for key_ in  pre_filt_kmers:
        if key_ in dict_back:
            del dict_foregr[key_]
