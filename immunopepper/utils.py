"""Contain functions to help compute, to preprocess"""
import itertools
import scipy as sp
import numpy as np

import sys
import bisect
import logging

from .constant import NOT_EXIST
from .immuno_nametuple import Coord
from .immuno_nametuple import Flag
from .immuno_nametuple import Idx
from .immuno_nametuple import OutputBackground
from .immuno_nametuple import Peptide
from .immuno_nametuple import ReadingFrameTuple


def to_adj_list(adj_matrix):
    """
    Converts a binary adjacency matrix to a list of directed edges
    """
    adj_list = []
    assert (adj_matrix.shape[0] == adj_matrix.shape[1])
    for idx in range(adj_matrix.shape[0]):
        for jdx in range(idx, adj_matrix.shape[0]):
            if adj_matrix[idx, jdx] == 1:
                adj_list.append([idx, jdx])
    return adj_list


def get_successor_list(adj_matrix, vertex_map, read_strand):
    """ Returns a list of successors by vertex, sensitive to the read strand

    Parameters
    ----------
    adj_matrix: 2D array with size [len(edges),len(edges)]. '1' represents the two edges
        are connected.
    vertex_map: 2D array with size [2, len(vertexs)]. Vertex positions
    read_strand: The gene strand. '+' or '-'.

    Returns
    -------
    succ_list: List[List[int]]. The connected vertex list for each vertex

    """
    succ_list = []
    assert (adj_matrix.shape[0] == adj_matrix.shape[1])
    for idx in range(adj_matrix.shape[0]):
        succ_list.append([])
        for jdx in range(adj_matrix.shape[0]):
            if adj_matrix[idx, jdx] == 1 and \
               ((read_strand == "+" and vertex_map[0][idx] <= vertex_map[0][jdx]) or \
                (read_strand == "-" and vertex_map[1][idx] >= vertex_map[1][jdx])):
                    succ_list[idx].append(jdx)
    return succ_list


def find_overlapping_cds_simple(v_start, v_stop, cds_begins, strand):
    """
    Find overlapping CDS within an exon given a list of CDS starts
    """
    # cds_start = cds_begin[0]
    if strand == '+':
        return list(filter(lambda x: x[0] >= v_start and x[0] < v_stop, cds_begins))
    else:
        return list(filter(lambda x: x[0] > v_start and x[0] <= v_stop, cds_begins))


def leq_strand(coord1, coord2, strand):
    if strand == "+":
        return coord1 <= coord2
    else:
        return coord1 >= coord2


def complementary_seq(dna_seq):
    """ Yields the complementary DNA sequence
    Only convert the character in comp_dict.
    Otherwise remain the same.
    """
    comp_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([comp_dict[nuc] if nuc in comp_dict else nuc for nuc in dna_seq])


def encode_chromosome(in_num):
    """
    Encodes chromosome to same cn
    """
    convert_dict = {23: "X", 24: "Y", 25: "MT"}
    return convert_dict[in_num] if in_num in convert_dict else str(in_num)


def translate_dna_to_peptide(dna_str):
    """ Translate a DNA sequence encoding a peptide to amino-acid sequence via RNA.

    If 'N' is included in input dna, 'X' will be outputted since 'N' represents
    uncertainty. Also will output a flag indicating if has stop codon.

    Parameters
    ----------
    dna_str: str or List(str). dna string to be translated.

    Returns
    -------
    aa_str: translated peptide
    has_stop_codon: Indicator for showing if the input dna contains stop codon

    """
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'
    }
    dna_str = dna_str.upper()
    has_stop_codon = False
    aa_str = []
    for idx in range(0, len(dna_str), 3):
        codon = dna_str[idx:idx + 3]
        if len(codon) < 3:
            break
        if 'N' in codon:
            aa_str.append('X')
        else:
            if codontable[codon] == '_':
                has_stop_codon = True
                return ''.join(aa_str), has_stop_codon
            else:
                aa_str.append(codontable[codon])

    return ''.join(aa_str), has_stop_codon


def get_sub_mut_dna(background_seq, coord, variant_comb, somatic_mutation_sub_dict, strand):
    """ Get the mutated dna sub-sequence according to mutation specified by the variant_comb.

    Parameters
    ----------
    background_seq: List(str). backgound sequence.
    start_v1: int. start position of first vertex.
    stop_v1: int. stop position of first vertex.
    start_v2: int. start position of second vertex.
    stop_v2: int. stop position of second vertex.
    variant_comb: List(int). List of variant position. Like ['38', '43']
    somatic_mutation_sub_dict: Dict. variant position -> variant details.
    strand: gene strand

    Returns
    -------
    final_dna: str. dna when applied somatic mutation.

    """
    def _get_variant_pos_offset(variant_pos, coord_pair_list, strand):
        offset = 0
        takes_effect = False
        for pair in coord_pair_list:
            if variant_pos in range(pair[0], pair[1]):
                if strand == '+':
                    offset += variant_pos - pair[0]
                else:
                    offset += pair[1] - variant_pos - 1
                takes_effect = True
                break
            else:
                offset = pair[1] - pair[0]

        return offset if takes_effect else NOT_EXIST

    real_coord = list(filter(lambda x: x != NOT_EXIST and x != None, coord))
    assert len(real_coord) % 2 == 0
    coord_pair_list = list(zip(real_coord[::2], real_coord[1::2]))

    if strand == '+':
        sub_dna = ''.join([background_seq[pair[0]:pair[1]] for pair in coord_pair_list])
    else:
        sub_dna = ''.join([background_seq[pair[0]:pair[1]][::-1] for pair in coord_pair_list])

    if variant_comb == NOT_EXIST : # no mutation exist
        return sub_dna

    relative_variant_pos = [_get_variant_pos_offset(variant_ipos, coord_pair_list, strand) for variant_ipos in variant_comb]
    for i,variant_ipos in enumerate(variant_comb):
        mut_base = somatic_mutation_sub_dict[variant_ipos]['mut_base']
        ref_base = somatic_mutation_sub_dict[variant_ipos]['ref_base']
        pos = relative_variant_pos[i]
        if pos != NOT_EXIST:
            sub_dna = sub_dna[:pos] + mut_base + sub_dna[pos+1:]
    return sub_dna


def cross_peptide_result(read_frame, strand, variant_comb, somatic_mutation_sub_dict, ref_mut_seq, peptide_accept_coord):
    """ Get translated peptide from the given exon pairs.

    Parameters
    ----------
    read_frame: NamedTuple. (read_start_codon, read_stop_codon, emitting_frame)
    strand: str. '+' or '-'
    variant_comb: List(int).
    somatic_mutation_sub_dict: Dict. variant position -> variant details.
    ref_mut_seq: Dict.['ref', 'background'] -> List(str)
    peptide_accept_coord: The start and end position of next vertex. Positions of the first vertex
        are already given in read_frame.

    Returns
    -------
    peptide: NamedTuple Peptide. has attribute ['ref', 'mut']. contain the output peptide
        translated from reference sequence and mutated sequence.
    coord: NamedTuple Coord. has attribute ['start_v1', 'stop_v1', 'start_v2', 'stop_v2']
        contains the true four position of exon pairs (after considering read framee)
        that outputs the peptide.
    flag: NamedTuple Flag. has attribute ['has_stop', 'is_isolated']
    next_reading_frame: Tuple. The reading frame to be propogated to the next vertex.

    """
    cds_left_modi, cds_right_modi, emitting_frame = read_frame[0], read_frame[1], read_frame[2]
    next_emitting_frame = (peptide_accept_coord[1] - peptide_accept_coord[0] + emitting_frame) % 3
    start_v1 = cds_left_modi
    stop_v1 = cds_right_modi

    #                                       |next_start_v1  |
    # |      v1           | |    v2                         |
    # -----[emitting_frame] [accepting_frame]-------
    # emitting_frame + accepting_frame = 3
    accepting_frame = (3 - emitting_frame) % 3

    if somatic_mutation_sub_dict:  # exist maf dictionary, so we use germline mutation-applied seq as the background seq
        ref_seq = ref_mut_seq['background']
    else:
        ref_seq = ref_mut_seq['ref']
    mut_seq = ref_mut_seq['background']
    # python is 0-based while gene annotation file(.gtf) is one based
    # so we need to do a little modification
    if strand == "+":
        start_v2 = peptide_accept_coord[0]
        stop_v2 = max(start_v2, peptide_accept_coord[1] - next_emitting_frame)
        coord = Coord(start_v1, stop_v1, start_v2, stop_v2)
        peptide_dna_str_mut = get_sub_mut_dna(mut_seq, coord, variant_comb, somatic_mutation_sub_dict, strand)
        peptide_dna_str_ref = ref_seq[start_v1:stop_v1] + ref_seq[start_v2:stop_v2]
        next_start_v1 = min(start_v2 + accepting_frame, peptide_accept_coord[1])
        next_stop_v1 = peptide_accept_coord[1]
    else:  # strand == "-"
        stop_v2 = peptide_accept_coord[1]
        start_v2 = min(stop_v2, peptide_accept_coord[0] + next_emitting_frame)
        coord = Coord(start_v1, stop_v1, start_v2, stop_v2)
        peptide_dna_str_mut = complementary_seq(get_sub_mut_dna(mut_seq, coord, variant_comb, somatic_mutation_sub_dict, strand))
        peptide_dna_str_ref = complementary_seq(ref_seq[start_v1:stop_v1][::-1] + ref_seq[start_v2:stop_v2][::-1])
        next_start_v1 = peptide_accept_coord[0]
        next_stop_v1 = max(stop_v2 - accepting_frame, peptide_accept_coord[0])

    next_reading_frame = ReadingFrameTuple(next_start_v1, next_stop_v1, next_emitting_frame)
    assert (len(peptide_dna_str_mut) == len(peptide_dna_str_ref))
    # if len(peptide_dna_str_mut) % 3 != 0:
    #     print("Applied mutations have changed the length of the DNA fragment - no longer divisible by 3")
    peptide_mut, mut_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_mut)
    peptide_ref, ref_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_ref)

    # if the stop codon appears before translating the second exon, mark 'single'
    is_isolated = False
    if len(peptide_mut)*3 <= abs(stop_v1 - start_v1) + 1:
        is_isolated = True
        jpos = 0.0
    else:
        jpos = float(stop_v1 - start_v1) / 3.0
    peptide = Peptide(peptide_mut, peptide_ref)
    coord = Coord(start_v1, stop_v1, start_v2, stop_v2)
    flag = Flag(mut_has_stop_codon, is_isolated)
    return peptide, coord, flag, next_reading_frame


def isolated_peptide_result(read_frame, strand, variant_comb, somatic_mutation_sub_dict, ref_mut_seq):
    """ Deal with translating isolated cds, almost the same as cross_peptide_result

    Parameters
    ----------
    read_frame: Tuple. (read_start_codon, read_stop_codon, emitting_frame)
    strand: str. '+' or '-'
    variant_comb: List(int).
    somatic_mutation_sub_dict: Dict. variant position -> variant details.
    ref_mut_seq: Dict.['ref', 'background'] -> List(str)

    Returns
    -------
    peptide: NamedTuple. has attribute ['ref', 'mut']. contain the output peptide
        translated from reference sequence and mutated sequence.
    coord: NamedTuple. has attribute ['start_v1', 'stop_v1', 'start_v2', 'stop_v2']
        contains the true two position of exon pairs (after considering read framee)
        that outputs the peptide. 'start_v2', 'stop_v2' is set to be NOT_EXIST.
    flag: NamedTuple. has attribute ['has_stop', 'is_isolated']

    """

    start_v1, stop_v1, emitting_frame = read_frame.cds_left_modi,read_frame.cds_right_modi,read_frame.read_phase
    start_v2 = NOT_EXIST
    stop_v2 = NOT_EXIST

    if somatic_mutation_sub_dict:  # exist maf dictionary, so we use germline mutation-applied seq as the background seq
        ref_seq = ref_mut_seq['background']
    else:
        ref_seq = ref_mut_seq['ref']
    mut_seq = ref_mut_seq['background']

    if strand == '+':
        coord = Coord(start_v1, stop_v1, start_v2, stop_v2)
        peptide_dna_str_mut = get_sub_mut_dna(mut_seq, coord, variant_comb, somatic_mutation_sub_dict, strand)
        peptide_dna_str_ref = ref_seq[start_v1:stop_v1]
    else:
        coord = Coord(start_v1, stop_v1, start_v2, stop_v2)
        peptide_dna_str_mut = complementary_seq(get_sub_mut_dna(mut_seq, coord, variant_comb, somatic_mutation_sub_dict, strand))
        peptide_dna_str_ref = complementary_seq(ref_seq[start_v1:stop_v1][::-1])

    peptide_mut, mut_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_mut)
    peptide_ref, ref_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_ref)

    is_isolated = True

    peptide = Peptide(peptide_mut,peptide_ref)
    coord = Coord(start_v1, stop_v1, start_v2, stop_v2)
    flag = Flag(mut_has_stop_codon,is_isolated)

    return peptide, coord, flag


def get_peptide_result(simple_meta_data, strand, variant_comb, somatic_mutation_sub_dict, ref_mut_seq):
    if somatic_mutation_sub_dict:  # exist maf dictionary, so we use germline mutation-applied seq as the background seq
        ref_seq = ref_mut_seq['background']
    else:
        ref_seq = ref_mut_seq['ref']
    mut_seq = ref_mut_seq['background']
    modi_coord = simple_meta_data.modified_exons_coord
    if strand == "+":
        peptide_dna_str_mut = get_sub_mut_dna(mut_seq,modi_coord, variant_comb, somatic_mutation_sub_dict, strand)
        peptide_dna_str_ref = get_sub_mut_dna(ref_seq,modi_coord, NOT_EXIST, somatic_mutation_sub_dict, strand)
    else:  # strand == "-"
        peptide_dna_str_mut = complementary_seq(get_sub_mut_dna(mut_seq, modi_coord, variant_comb, somatic_mutation_sub_dict, strand))
        peptide_dna_str_ref = complementary_seq(get_sub_mut_dna(ref_seq, modi_coord, NOT_EXIST, somatic_mutation_sub_dict, strand))
    peptide_mut, mut_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_mut)
    peptide_ref, ref_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_ref)

    # if the stop codon appears before translating the second exon, mark 'single'
    stop_v1 = modi_coord.stop_v1
    start_v1 = modi_coord.start_v1
    start_v2 = modi_coord.start_v2
    if start_v2 == NOT_EXIST or len(peptide_mut)*3 <= abs(stop_v1 - start_v1) + 1:
        is_isolated = True
    else:
        is_isolated = False
    peptide = Peptide(peptide_mut,peptide_ref)
    flag = Flag(mut_has_stop_codon,is_isolated)
    return peptide,flag


def get_size_factor(strains, lib_file_path):
    """ Get the expression adjustment parameter for certain samples.

    Parameters
    ----------
    samples: List[str]. List of samples.
    lib_file_path: str. libsize file path

    Returns
    -------
    sf: 1darray. size factor.

    """
    libs = sp.loadtxt(lib_file_path, dtype='str', skiprows=1, delimiter='\t')
    a, b = sp.where(strains[:, sp.newaxis] == libs[:, 0])
    assert sp.all(libs[b, 0] == strains)
    libs = libs[b, :]
    med = sp.median(libs[:, 1].astype('float'))
    sf = med / libs[:, 1].astype('float')
    return sf


def get_all_comb(array, r=None):
    """ Get all combinations of items in the given array
    Specifically used for generating variant combination

    Parameters
    ----------
    array: 1D array. input array
    r: int. The number of items in a combination

    Returns
    -------
    result: List(Tuple). List of combination

    """
    if r is None:
        r = len(array)
    result = []
    for i in range(1,r+1):
        result.extend(list(itertools.combinations(array,i)))
    return result


def is_mutation_within_exon_pair(pos_start, pos_end, mutation_pos):
    #variant_pos = [pos for pos in mutation_pos]
    #variant_pos_candi = [ipos for ipos in variant_pos if ipos > pos_start and ipos < pos_end]
    return list(filter(lambda x: x > pos_start and x < pos_end, mutation_pos))


def is_isolated_cds(gene, idx):
    """
    Indicate if a peptide is translated by an isolated cds

    Currently, it is a simple version. If the exon where the cds is
    has no connected exon, it is isolated. (the method is not true in some cases)
    Parameters
    ----------
    gene: gene object
    idx: exon id

    """

    if len(gene.vertex_succ_list[idx]) > 0:
        return False

    return sp.sum(gene.splicegraph.edges[:, idx]) == 0


def is_in_junction_list(v1,v2,strand,junction_list):
    """Check if the intron is in concerned junction list"""
    if junction_list is not None:
        return int(':'.join([str(v1[1]),  str(v2[0]), strand]) in junction_list)
    else:
        return NOT_EXIST


def get_exon_expr(gene,vstart,vstop,Segments,Idx):
    """ Split the exon into segments and get the corresponding counts.

    Parameters
    ----------
    gene: Object. Generated by SplAdder.
    vstart: int. Start position of given exon.
    vstop: int. Stop position of given exon.
    Segments: Namedtuple, store segment expression information from count.hdf5.
        has attribute ['expr', 'lookup_table'].
    Idx: Namedtuple, has attribute idx.gene and idx.sample

    Returns
    -------
    expr_list: List[Tuple(int,float)] (int, float) represents the length of segment
        and the expression count of that segment.

    """
    # Todo: deal with absense of count file
    if vstart == NOT_EXIST or vstop == NOT_EXIST:  # isolated exon case
        expr_list = []
        return expr_list
    if Segments is None or Idx.sample is None:
        expr_list = [NOT_EXIST]
        return expr_list
    count_segments = Segments.lookup_table[gene.name]
    segments = gene.segmentgraph.segments
    seg_count = Segments.expr[count_segments,Idx.sample]
    sv1_id = bisect.bisect(segments[0], vstart)-1
    sv2_id = bisect.bisect(segments[0], vstop)-1
    expr_list = []
    if sv1_id == sv2_id:
        expr_list.append((vstop - vstart, seg_count[sv1_id]))
    else:
        expr_list.append((segments[1, sv1_id] - vstart, seg_count[sv1_id]))
        for i in range(sv1_id + 1, sv2_id):
            expr_list.append((segments[1, i] - segments[0, i], seg_count[i]))
        expr_list.append((vstop-segments[0, sv2_id], seg_count[sv2_id]))
        if gene.strand == '-': # need to reverse epression list to match the order of translation
            expr_list = expr_list[::-1]
    return expr_list


def get_segment_expr(gene, coord, Segments, Idx):
    """ Get the segment expression for one exon-pair.
    Apply 'get_exon_expr' for each exon and concatenate them.

    Parameters
    ----------
    gene: Object. Generated by SplAdder.
    coord: NamedTuple. has attribute ['start_v1', 'stop_v1', 'start_v2', 'stop_v2']
        contains the true four position of exon pairs (after considering read framee)
        that outputs the peptide.
    Segments: Namedtuple, store segment expression information from count.hdf5.
        has attribute ['expr', 'lookup_table'].
    Idx: Namedtuple, has attribute idx.gene and idx.sample

    Returns
    -------
    mean_expr: float. The average expression counts for the given exon-pair
    expr1: List[Tuple(int,float)] (int, float) represents the length of segment
        and the expression count of that segment.

    """
    expr_list1 = get_exon_expr(gene,coord.start_v1,coord.stop_v1,Segments,Idx)
    expr_list2 = get_exon_expr(gene,coord.start_v2,coord.stop_v2,Segments,Idx)
    expr_list1.extend(expr_list2)
    if coord.start_v3 is not None:
        expr_list3 = get_exon_expr(gene,coord.start_v3,coord.stop_v3,Segments,Idx)
        expr_list1.extend(expr_list3)
    expr_sum = 0
    seg_len = 0
    for item in expr_list1:
        length = item[0]
        expr = item[1]
        expr_sum += length*expr
        seg_len += length
    mean_expr = int(expr_sum/seg_len) if seg_len > 0 else 0
    return mean_expr,expr_list1


def get_total_gene_expr(gene, Segments, Idx):
    """ get total reads count for the given sample and the given gene
    actually total_expr = reads_length*total_reads_counts
    """
    if Segments is None or Idx.sample is None:
        return NOT_EXIST
    count_segments = Segments.lookup_table[gene.name]
    seg_len = gene.segmentgraph.segments[1]-gene.segmentgraph.segments[0]
    seg_expr = Segments.expr[count_segments,Idx.sample]
    total_expr = np.sum(seg_len*seg_expr)
    return total_expr


def get_idx(sample_idx_table, sample, gene_idx):
    """ Create a aggregated Index with nametuple idx
    Combine the gene_idx, sample_idx

    Parameters
    ----------
    sample_idx_table: Dict. str -> int. Mapping sample to idx.
    sample: str. sample name
    gene_idx: int. Gene index, mainly for formatting the output.

    """
    if sample_idx_table is not None:
        if sample in sample_idx_table:
            sample_idx = sample_idx_table[sample]
        else:
            sample_idx = None
            logging.warning("utils.py: The sample {} is not in the count file. Program proceeds without outputting expression data.".format(sample))
    else:
        sample_idx = None
    idx = Idx(gene_idx,sample_idx)
    return idx


def create_libsize(expr_distr_dict,output_fp,debug=False):
    """ create library_size text file.

    Calculate the 75% expression and sum of expression for each sample
    and write into output_fp.

    Parameters
    ----------
    expr_distr_dict: Dict. str -> List(float). Mapping sample to the expression of all exon pairs
    output_fp: file pointer. library_size text
    debug: Bool. In debug mode, return the libsize_count dictionary.
    """
    # filter the dict
    libsize_count = {}
    for sample,expr_list in expr_distr_dict.items():
        if np.array(expr_list).dtype in  [np.float,np.int]:
            libsize_count[sample] = (np.percentile(expr_list,75),np.sum(expr_list))
    #libsize_count = {sample:(np.percentile(expr_list,75),np.sum(expr_list)) for sample,expr_list in list(expr_distr_dict.items())}
    if debug:
        return libsize_count
    with open(output_fp,'w') as f:
        f.write('\t'.join(['sample','libsize_75percent','libsize_total_count'])+'\n')
        for sample,count_tuple in list(libsize_count.items()):
            line = '\t'.join([sample,str(round(count_tuple[0],1)),str(int(count_tuple[1]))])+'\n'
            f.write(line)

def get_concat_peptide(front_coord_pair, back_coord_pair,front_peptide, back_peptide, strand,k=None):
    """
    Get the concatenated peptide from possible match peptide.
    Parameters
    ----------
    front_coord_pair: str. coordinate string for the front exon pair.
    back_coord_pair: str. coordinate string for the back exon pair.
    front_peptide: str. peptide translated from the front exon pair.
    back_peptide: str. peptide translated from the back exon pair.
    strand: str. '+' or '-'
    k: k for k-mer.

    Returns
    -------
    new_peptide: str. concatenated peptide. If none is found, return empty string

    Examples
    --------
    front_pep: 'MKTT', back_pep: 'TTAC', concatenated_pep: 'MKTTAC'

    """
    def get_longest_match_position(front_str,back_str,L=None):
        if L is None:
            L = min(len(front_str),len(back_str))
        for i in reversed(list(range(1,L+1))):
            if front_str[-i:] == back_str[:i]:
                return i
        return None
    if strand == '+':
        front_coord = front_coord_pair.stop_v2
        back_coord = back_coord_pair.start_v1
    else:
        front_coord = front_coord_pair.start_v2
        back_coord = back_coord_pair.stop_v1
    if abs(front_coord-back_coord) % 3 == 0:
        if front_coord == back_coord:  # no intersection and we concatenate them directly
            new_peptide = front_peptide + back_peptide
        else:
            pep_common_num = get_longest_match_position(front_peptide,back_peptide,L=k)
            if pep_common_num is None:
                new_peptide = ''
            else:
                new_peptide = front_peptide + back_peptide[pep_common_num:]
        return new_peptide
    else:
        return ''

def get_concat_junction_peptide(gene, output_peptide_list, output_metadata_list, Segments, Idx, k):
    '''
    Find all the match peptide and concatenate them, output the new peptide and its expression list

    Parameters
    ----------
    gene: SplAdder object.
    output_peptide_list: List[namedtuple]. Contain all the possible output peptide in the given splicegraph.
    output_metadata_list: List[namedtuple]. Contain the correpsonding medata data for each output peptide.
    Idx: Namedtuple Idx, has attribute idx.gene and idx.sample
    Segments: Namedtuple Segments, store segment expression information from count.hdf5.
           has attribute ['expr', 'lookup_table'].
    k: k for k-mer. Positive k is required so that we can find key vertex

    Returns
    -------
    concat_peptide_list: List[str]. Contain all the possible concatenated peptide in the given splicegraph.
    concat_expr_list:List[List(Tuple(int,float))]. Contain the segment expression data for each concatenated peptide.

    '''
    def get_concat_expr_list(front_coord_pair,back_vertex_pair):
        new_expr_list = []
        # first vertex

        expr_list = get_exon_expr(gene, front_coord_pair.start_v1, front_coord_pair.stop_v1, Segments, Idx)
        new_expr_list.extend(expr_list)
        # second vetex
        back_start,back_end = back_vertex_pair.split(',')
        second_cds_id = int(back_start)
        expr_list = get_exon_expr(gene, vertices[0, second_cds_id], vertices[1, second_cds_id], Segments, Idx)
        new_expr_list.extend(expr_list)
        # third vertex
        if back_end == '.':
            expr_list = [(0,0)]
        else:
            third_cds_id = int(back_end)
            expr_list = get_exon_expr(gene, vertices[0, third_cds_id], vertices[1, third_cds_id], Segments, Idx)
        new_expr_list.extend(expr_list)
        return new_expr_list
    vertices = gene.splicegraph.vertices
    vertex_len = vertices[1,:]-vertices[0,:]
    key_id_list = np.where(vertex_len < (k+1)*3)[0]

    vertex_id_pair_list = [metadata.vertex_idx for metadata in output_metadata_list]
    coord_pair_list = [metadata.exons_coor for metadata in output_metadata_list]
    concat_peptide_list = []
    concat_expr_lists = []
    for key_id in key_id_list:
        front_id_list =[i for i, vert_pair in enumerate(vertex_id_pair_list) if vert_pair[1] == str(key_id)]
        back_id_list =[i for i, vert_pair in enumerate(vertex_id_pair_list) if vert_pair[0] == str(key_id)]
        for front_id in front_id_list:
            for back_id in back_id_list:
                triple_v = '_'.join([vertex_id_pair_list[front_id][0],vertex_id_pair_list[front_id][1],vertex_id_pair_list[back_id][1]])
                back_peptide = output_peptide_list[back_id].peptide
                front_peptide = output_peptide_list[front_id].peptide
                back_coord_pair = coord_pair_list[back_id]
                front_coord_pair = coord_pair_list[front_id]
                if len(back_peptide) > 0 and len(front_peptide) > 0 and back_coord_pair[-1] != NOT_EXIST: # filter out those empty string
                    concat_peptide = get_concat_peptide(front_coord_pair,back_coord_pair,front_peptide, back_peptide,gene.strand,k)
                    if len(concat_peptide) > 0 : # calculate expr list
                        concat_expr_list = get_concat_expr_list(front_coord_pair,vertex_id_pair_list[back_id])
                        concat_expr_lists.append(concat_expr_list)
                        concat_peptide = OutputBackground(id='Gene'+str(Idx.gene)+'_'+triple_v,
                                                           peptide=concat_peptide)
                        concat_peptide_list.append(concat_peptide)
    return concat_peptide_list,concat_expr_lists

def convert_namedtuple_to_str(_namedtuple,field_list):
    def convert_list_to_str(_list):
        remove_none_list = filter(lambda x:x is not None, _list)
        return ';'.join([str(_item) for _item in remove_none_list])
    line = ''
    for field in field_list:
        if field == 'new_line':
            line = line.strip()+'\n'
            continue
        # should first check if field in namedtuple
        item = getattr(_namedtuple, field)
        if isinstance(item,(list,tuple)):
            line += convert_list_to_str(item)+'\t'
        else:
            line += str(item)+'\t'
    return line[:-1] # remove the last '\t'

def write_namedtuple_list(fp, namedtuple_list, field_list):
    """ Write namedtuple_list to the given file pointer"""
    fp.writelines(convert_namedtuple_to_str(_namedtuple,field_list)+'\n' for _namedtuple in namedtuple_list)

def write_list(fp, _list):
    fp.writelines([l+'\n' for l in _list])

def write_gene_expr(fp, gene_expr_tuple_list):
    header_line = 'gene\texpr\n'
    fp.write(header_line)
    gene_expr_str_list = [ gene_expr_tuple[0]+'\t'+str(gene_expr_tuple[1]) for gene_expr_tuple in gene_expr_tuple_list]
    write_list(fp,gene_expr_str_list)


def check_chr_consistence(ann_chr_set,mutation,graph_data):
    germline_chr_set = set()
    somatic_chr_set = set()
    mode = mutation.mode
    if mutation.germline_mutation_dict:
        germline_chr_set = set([item[1] for item in mutation.germline_mutation_dict.keys()])
    if mutation.somatic_mutation_dict:
        somatic_chr_set = set([item[1] for item in mutation.somatic_mutation_dict.keys()])
    whole_mut_set = somatic_chr_set.union(germline_chr_set)
    common_chr = whole_mut_set.intersection(ann_chr_set)

    if mode != 'ref' and len(common_chr) == 0:
        logging.error("Mutation file has different chromosome naming from annotation file, please check")
        sys.exit(0)
    if len(graph_data) > 0:
        gene_chr_set = set([gene.chr for gene in graph_data])
        new_chr_set = gene_chr_set.difference(ann_chr_set)
        if len(new_chr_set) > 0:
            logging.error("Gene object has different chromosome naming from annotation file, please check")
            sys.exit(0)

def codeUTF8(s):
    return s.encode('utf-8')

def decodeUTF8(s):
    if not hasattr(s, 'decode'):
        return s
    return s.decode('utf-8')
