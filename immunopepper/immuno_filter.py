"""Contain the function used for filtering"""
import sys
import re

import scipy as sp

from .utils import complementary_seq,translate_dna_to_peptide,get_exon_expr
from .immuno_nametuple import OutputBackground
from .constant import NOT_EXIST
from functools import reduce

def junction_is_annotated(gene, gene_to_transcript_table, transcript_to_cds_table):
    """ Indicate whether exon pair also appears in transcript given by .gtf file

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

    junction_flag = sp.zeros(edges.shape, dtype='bool')
    for ts in ts_list:
        if not ts in transcript_to_cds_table:
            continue

        curr_ts = transcript_to_cds_table[ts]
        if len(curr_ts) < 2:
            continue
        curr_ts = sp.array(curr_ts)[:, [0, 1]]
        transcript_junctions = [':'.join(x) for x in
                                curr_ts.ravel()[1:-1].reshape(curr_ts.shape[0] - 1, 2).astype('str')]

        for x in sp.array(sp.where(sp.triu(edges))).T:
            if '%i:%i' % (vertices[1, x[0]], vertices[0, x[1]]) in transcript_junctions:
                junction_flag[x[0], x[1]] = 1
                junction_flag[x[1], x[0]] = 1

    return junction_flag

def get_junction_anno_flag(junction_flag, vertex_id_tuple):
    """
    Output is_junction_annotated flag
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
        return NOT_EXIST
    junction_list = [(vertex_id_tuple[i], vertex_id_tuple[i+1]) for i in range(len(vertex_id_tuple)-1)]
    if len(junction_list) > 1:
        #return reduce(lambda junction1,junction2: junction_flag[junction1[0],junction1[1]]+junction_flag[junction2[0],junction2[1]],junction_list)
        return list(map(lambda junction: int(junction_flag[junction[0],junction[1]]),junction_list))
    else:
        return int(junction_flag[vertex_id_tuple[0],vertex_id_tuple[1]])

def get_full_peptide(gene, seq, cds_list, Segments, Idx, mode):
    """
    Output translated peptide and segment expression list given cds_list
    Parameters
    ----------
    gene: Object, created by SplAdder.
    seq: str. Gene sequence.
    cds_list: List[Tuple(v_start,v_stop,reading_frame)]
    Segments: Namedtuple Segments, store segment expression information from count.hdf5.
           has attribute ['expr', 'lookup_table'].
    Idx: Namedtuple Idx, has attribute idx.gene and idx.sample
    mode: [temporal argument]. Due to the different meaning of cds_tuple in gene.splicegraph.reading_frame
        and that in gene_to_cds dict (the formal v_start and v_stop has already considered reading_frame and
        do not need additional modified), we need to deal with them differently. So for now there are two modes
        'full' and 'background', indicating full-kmer and background. Will remove it in the future version.

    Returns
    -------
    cds_expr_list: List[Tuple(segment_length,segment_expression)]
    cds_string: str. Concatenated sequence string according to cds_list
    cds_peptide: str. Translated peptide string according to cds_list

    """
    if gene.strand.strip() == "-" and mode=='back':
        cds_list = cds_list[::-1]

    cds_string = ""
    first_cds = True
    cds_expr_list = []
    # Append transcribed CDS regions to the output
    for coord_left, coord_right, frameshift in cds_list:

        # Apply initial frameshift on the first CDS of the transcript
        if first_cds and mode != 'full':
            if gene.strand.strip() == "+":
                coord_left += frameshift
            else:
                coord_right -= frameshift
            first_cds = False
        cds_expr = get_exon_expr(gene, coord_left, coord_right, Segments, Idx)
        cds_expr_list.extend(cds_expr)
        nuc_seq = seq[coord_left:coord_right]

        # Accumulate new DNA sequence...
        if gene.strand.strip() == "+":
            cds_string += nuc_seq
        elif gene.strand.strip() == "-":
            cds_string += complementary_seq(nuc_seq[::-1])
        else:
            print("ERROR: Invalid strand. Got %s but expect + or -" % gene.strand.strip())
            sys.exit(1)
    cds_peptide, is_stop_flag = translate_dna_to_peptide(cds_string)
    return cds_expr_list, cds_string, cds_peptide


def peptide_match(background_peptide_list, peptide):
    """ Find if the translated exon-pair peptide also appear in the background peptide translated from annotation file.

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
        transcript_id = background_peptide.id
        if not re.search(peptide, ref_peptide) is None:
            match_ts_list.append(transcript_id)
    return match_ts_list


def get_exon_dict(metadata_list,strand):
    """ Get an auxiliary dictionary used for filtering

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
        coord = metadata.exons_coor
        idx = metadata.output_id
        read_frame = metadata.read_frame.read_phase
        if strand == '+':
            key = (read_frame,coord.stop_v1,coord.start_v2)
            if key in exon_dict:
                exon_dict[key].append((idx, coord.start_v1, coord.stop_v2))
            else:
                exon_dict[key] = [(idx, coord.start_v1, coord.stop_v2)]
        else:
            key = (read_frame, coord.stop_v2, coord.start_v1)
            if key in exon_dict:
                exon_dict[key].append((idx, coord.start_v2, coord.stop_v1))
            else:
                exon_dict[key] = [(idx, coord.start_v2, coord.stop_v1)]
    return exon_dict


def get_remove_id(metadata_dict):
    # TODO: [Warning!] if two output lines have identical
    # coordinates and readframe
    # both of them will be removed. Need to fix.
    remove_id_list = []
    for exon_pair_list in list(metadata_dict.values()):
        L = len(exon_pair_list)
        if L < 2:
            continue
        for i, exon_pair in enumerate(exon_pair_list):
            for j in range(L):
                i_pos1 = exon_pair[1]
                i_pos2 = exon_pair[2]
                j_pos1 = exon_pair_list[j][1]
                j_pos2 = exon_pair_list[j][2]
                if j!=i and i_pos1 >= j_pos1 and i_pos2 <= j_pos2:
                    remove_id_list.append(exon_pair[0])
                    break
    return remove_id_list

def get_filtered_metadata_list(metadata_list,strand):
    exon_dict = get_exon_dict(metadata_list,strand)
    remove_id_list = get_remove_id(exon_dict)
    filtered_meta_list = list(filter(lambda metadata: metadata.output_id not in remove_id_list,metadata_list))
    return filtered_meta_list


def get_filtered_output_list(metadata_list,peptide_list,expr_lists):
    """ Get all redundant output id.

    Parameters
    ----------
    metadata_list: List[Output_metadata].
    peptide_list: List[Output_junc_peptide]
    expr_lists:: List[List(tuple)]

    Returns
    -------
    remove_id_list: List[str]. The list of id to be removed.
    """
    exon_dict = get_exon_dict(metadata_list)
    remove_id_list = get_remove_id(exon_dict)
    filtered_meta_list = []
    filtered_peptide_list = []
    filtered_expr_lists = []
    for i,metadata in enumerate(metadata_list):
        idx = metadata.output_id
        if idx not in remove_id_list:
            filtered_meta_list.append(metadata)
            filtered_peptide_list.append(peptide_list[i])
            filtered_expr_lists.append(expr_lists[i])
    return filtered_meta_list, filtered_peptide_list,filtered_expr_lists



