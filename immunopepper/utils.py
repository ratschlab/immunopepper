"""Contain functions to help compute, to preprocess"""
import itertools
import scipy as sp
import numpy as np
from anytree import Node, RenderTree

from collections import namedtuple
import bisect

from constant import NOT_EXIST

Peptide = namedtuple('Peptide', ['mut', 'ref'])
Coord = namedtuple('Coord', ['start_v1', 'stop_v1', 'start_v2', 'stop_v2'])
Flag = namedtuple('Flag', ['has_stop', 'is_isolated'])
Idx = namedtuple('Idx', ['gene', 'sample'])


def to_adj_list(adj_matrix):
    """
    Converts a binary adjacency matrix to a list of directed edges
    """
    adj_list = []
    assert (adj_matrix.shape[0] == adj_matrix.shape[1])
    for idx in range(adj_matrix.shape[0]):
        for jdx in range(adj_matrix.shape[0]):
            if adj_matrix[idx, jdx] == 1 and idx <= jdx:
                adj_list.append([idx, jdx])
    return adj_list


def to_adj_succ_list(adj_matrix, vertex_map, read_strand):
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
    return filter(lambda cds_begin: cds_begin[0] >= v_start and cds_begin[0] <= v_stop, cds_begins)


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
    comp_dict_keys = comp_dict.keys()
    return "".join(map(lambda nuc: comp_dict[nuc] if nuc in comp_dict_keys else nuc, dna_seq))


def encode_chromosome(in_num):
    """
    Encodes chromosome to same cn
    """
    convert_dict = {23: "X", 24: "Y", 25: "MT"}
    return convert_dict[in_num] if in_num in convert_dict.keys() else str(in_num)


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
        if "N" in codon:
            aa_str.append('X')
        else:
            if codontable[codon] == '_':
                has_stop_codon = True
                return ''.join(aa_str),has_stop_codon
            else:
                aa_str.append(codontable[codon])

    return "".join(aa_str),has_stop_codon


def get_sub_mut_dna(background_seq, start_v1, stop_v1, start_v2, stop_v2, variant_comb, mutation_sub_dic_maf, strand):
    """ Get the mutated dna sub-sequence according to mutation specified by the variant_comb.

    Parameters
    ----------
    background_seq: List(str). backgound sequence.
    start_v1: int. start position of first vertex.
    stop_v1: int. stop position of first vertex.
    start_v2: int. start position of second vertex.
    stop_v2: int. stop position of second vertex.
    variant_comb: List(int). List of variant position. Like ['38', '43']
    mutation_sub_dic_maf: Dict. variant position -> variant details.
    strand: gene strand

    Returns
    -------
    final_dna: str. dna when applied somatic mutation.

    """
    if start_v2 != NOT_EXIST:
        if strand == '-':
            sub_dna_list = list(background_seq[start_v1:stop_v1][::-1] + background_seq[start_v2:stop_v2][::-1])
        else:
            sub_dna_list = list(background_seq[start_v1:stop_v1] + background_seq[start_v2:stop_v2])
    else:
        if strand == '-':
            sub_dna_list = list(background_seq[start_v1:stop_v1][::-1])
        else:
            sub_dna_list = list(background_seq[start_v1:stop_v1])

    if variant_comb == NOT_EXIST : # no mutation exist
        return ''.join(sub_dna_list)
    for variant_ipos in variant_comb:
        mut_base = mutation_sub_dic_maf[variant_ipos]['mut_base']
        ref_base = mutation_sub_dic_maf[variant_ipos]['ref_base']
        # strand = mutation_sub_dic_maf[variant_ipos]['strand']

        # decide mutation happens in which exon vertice
        # it may falls out of the two ranges due to readframe shift
        if variant_ipos in range(start_v1, stop_v1):

            if strand == '-':
                pos = stop_v1-variant_ipos-1
                assert (sub_dna_list[pos] == ref_base)
                sub_dna_list[pos] = mut_base
            else:
                pos = variant_ipos - start_v1
                assert (sub_dna_list[pos] == ref_base)
                sub_dna_list[pos] = mut_base

        elif start_v2 != NOT_EXIST and variant_ipos in range(start_v2, stop_v2):
            if strand == '-':
                pos = stop_v2-variant_ipos+stop_v1-start_v1-1
                assert (sub_dna_list[pos] == ref_base)
                sub_dna_list[pos] = mut_base
            else:
                pos = variant_ipos-start_v2+stop_v1-start_v1
                assert (sub_dna_list[pos] == ref_base)
                sub_dna_list[pos] = mut_base

    final_dna = ''.join(sub_dna_list)
    return final_dna


def cross_peptide_result(read_frame, strand, variant_comb, mutation_sub_dic_maf, ref_mut_seq, peptide_accept_coord):
    """ Get translated peptide from the given exon pairs.

    Parameters
    ----------
    read_frame: Tuple. (read_start_codon, read_stop_codon, emitting_frame)
    strand: str. '+' or '-'
    variant_comb: List(int).
    mutation_sub_dic_maf: Dict. variant position -> variant details.
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
    cds_left_modi, cds_right_modi, emitting_frame = read_frame
    next_emitting_frame = (peptide_accept_coord[1] - peptide_accept_coord[0] + emitting_frame) % 3
    start_v1 = cds_left_modi
    stop_v1 = cds_right_modi

    #                                       |next_start_v1  |
    # |      v1           | |    v2                         |
    # -----[emitting_frame] [accepting_frame]-------
    # emitting_frame + accepting_frame = 3
    accepting_frame = (3 - emitting_frame) % 3

    if mutation_sub_dic_maf is None:
        ref_seq = ref_mut_seq['ref']
    else:
        ref_seq = ref_mut_seq['background']
    # python is 0-based while gene annotation file(.gtf) is one based
    # so we need to do a little modification
    if strand == "+":
        start_v2 = peptide_accept_coord[0]
        stop_v2 = peptide_accept_coord[1] - next_emitting_frame
        peptide_dna_str_mut = get_sub_mut_dna(ref_mut_seq['background'], start_v1, stop_v1, start_v2, stop_v2, variant_comb, mutation_sub_dic_maf, strand)
        peptide_dna_str_ref = ref_seq[start_v1:stop_v1] + ref_seq[start_v2:stop_v2]
        next_start_v1 = start_v2 + accepting_frame
        next_stop_v1 = peptide_accept_coord[1]
    else:  # strand == "-"
        start_v2 = peptide_accept_coord[0] + next_emitting_frame
        stop_v2 = peptide_accept_coord[1]
        mut_seq = get_sub_mut_dna(ref_mut_seq['background'], start_v1, stop_v1, start_v2, stop_v2, variant_comb, mutation_sub_dic_maf, strand)
        peptide_dna_str_mut = complementary_seq(mut_seq)
        peptide_dna_str_ref = complementary_seq(ref_seq[start_v1:stop_v1][::-1] + ref_seq[start_v2:stop_v2][::-1])
        next_start_v1 = peptide_accept_coord[0]
        next_stop_v1 = stop_v2 - accepting_frame

    next_reading_frame = (next_start_v1, next_stop_v1, next_emitting_frame)
    assert (len(peptide_dna_str_mut) == len(peptide_dna_str_ref))
    if len(peptide_dna_str_mut) % 3 != 0:
        print("Applied mutations have changed the length of the DNA fragment - no longer divisible by 3")
    peptide_mut, mut_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_mut)
    peptide_ref, ref_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_ref)

    # if the stop codon appears before translating the second exon, mark 'single'
    is_isolated = False
    if len(peptide_mut)*3 <= abs(stop_v1 - start_v1) + 1:
        is_isolated = True
        jpos = 0.0
    else:
        jpos = float(stop_v1 - start_v1) / 3.0
    peptide = Peptide(peptide_mut,peptide_ref)
    coord = Coord(start_v1,stop_v1,start_v2,stop_v2)
    flag = Flag(mut_has_stop_codon,is_isolated)
    return peptide, coord, flag, next_reading_frame


def isolated_peptide_result(read_frame, strand, variant_comb, mutation_sub_dic_maf, ref_mut_seq):
    """ Deal with translating isolated cds, almost the same as cross_peptide_result

    Parameters
    ----------
    read_frame: Tuple. (read_start_codon, read_stop_codon, emitting_frame)
    strand: str. '+' or '-'
    variant_comb: List(int).
    mutation_sub_dic_maf: Dict. variant position -> variant details.
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
    Peptide = namedtuple('Peptide',['mut','ref'])
    Coord = namedtuple('Coord',['start_v1','stop_v1','start_v2','stop_v2'])
    Flag = namedtuple('Flag', ['has_stop', 'is_isolated'])

    start_v1, stop_v1, emitting_frame = read_frame
    start_v2 = NOT_EXIST
    stop_v2 = NOT_EXIST

    if mutation_sub_dic_maf is None:  # no somatic mutation, the germline mutation will be the spotlight
        ref_seq = ref_mut_seq['ref']
    else:
        ref_seq = ref_mut_seq['background']  # we focus on somatic mutation in this case
    mut_seq = ref_mut_seq['background']

    if strand == '+':
        peptide_dna_str_mut = get_sub_mut_dna(mut_seq, start_v1, stop_v1, start_v2, stop_v2, variant_comb, mutation_sub_dic_maf, strand)
        peptide_dna_str_ref = ref_seq[start_v1:stop_v1]
    else:
        peptide_dna_str_mut = complementary_seq(get_sub_mut_dna(mut_seq, start_v1, stop_v1, start_v2, stop_v2, variant_comb, mutation_sub_dic_maf, strand))
        peptide_dna_str_ref = complementary_seq(ref_seq[start_v1:stop_v1][::-1])

    peptide_mut, mut_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_mut)
    peptide_ref, ref_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_ref)

    is_isolated = True

    peptide = Peptide(peptide_mut,peptide_ref)
    coord = Coord(start_v1,stop_v1,start_v2,stop_v2)
    flag = Flag(mut_has_stop_codon,is_isolated)

    return peptide, coord, flag


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
    """ Get all the combination of items in the given array
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


def is_mutation_within_exon_pair(pos_start,pos_end,mutation_pos):
    variant_pos = [pos for pos in mutation_pos]
    variant_pos_candi = [ipos for ipos in variant_pos if ipos > pos_start and ipos < pos_end]
    return variant_pos_candi


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
    return int(':'.join([str(v1[1]),  str(v2[0]), strand]) in junction_list)


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
    if vstart == NOT_EXIST or vstop == NOT_EXIST:  # isolated exon case
        expr_list = []
        return expr_list
    segments = gene.segmentgraph.segments
    sv1_id = bisect.bisect(segments[0], vstart)-1
    sv2_id = bisect.bisect(segments[0], vstop)-1
    expr_list = []
    if sv1_id == sv2_id:
        expr_list.append((vstop - vstart, Segments.expr[sv1_id, Idx.sample]))
    else:
        expr_list.append((segments[1, sv1_id] - vstart, Segments.expr[sv1_id,Idx.sample]))
        for i in range(sv1_id + 1, sv2_id):
            expr_list.append((segments[1, i] - segments[0, i], Segments.expr[i,Idx.sample]))
        expr_list.append((vstop-segments[0, sv2_id], Segments.expr[sv2_id, Idx.sample]))
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
    expr_sum = 0
    seg_len = 0
    for item in expr_list1:
        length = item[0]
        expr = item[1]
        expr_sum += length*expr
        seg_len += length
    mean_expr = int(expr_sum/seg_len)
    return mean_expr,expr_list1

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
        sample_idx = sample_idx_table[sample]
    else:
        sample_idx = None
    idx = Idx(gene_idx,sample_idx)
    return idx


def create_libsize(expr_distr_dict,output_fp):
    """ create library_size text file.

    Calculate the 75% expression and sum of expression for each sample
    and write into output_fp.

    Parameters
    ----------
    expr_distr_dict: Dict. str -> List(float). Mapping sample to the expression of all exon pairs
    output_fp: file pointer. library_size text

    """
    libsize_count = {sample:(np.percentile(expr_list,75),np.sum(expr_list)) for sample,expr_list in expr_distr_dict.items()}
    with open(output_fp,'w') as f:
        f.write('\t'.join(['sample','libsize_75percent','libsize_total_count'])+'\n')
        for sample,count_tuple in libsize_count.items():
            line = '\t'.join([sample,str(round(count_tuple[0],1)),str(int(count_tuple[1]))])+'\n'
            f.write(line)

def build_kmer_dict(kmer_file):
    f = open(kmer_file,'r')
    kmer_dict = {}
    for line in f:
        items = line.strip().split('\t')
        kmer = items[0]
        expr = float(items[1])
        if kmer not in kmer_dict:
            kmer_dict[kmer] = expr
        elif expr > kmer_dict[kmer]:
            kmer_dict[kmer] = expr
    return kmer_dict


def get_all_paths(gene):
    """
    Get all possible paths according to reading frames and success vertex list
    Parameters
    ----------
    gene: SpAddler gene object.

    Returns
    -------
    path_dict: Dict. Key is vertex id, value is the list of paths that end with that vertex.

    """
    def add_v_to_list3(_list,v):
        return [_ilist + [v] for _ilist in _list]
    reading_frame = gene.splicegraph.reading_frames
    succ_list = gene.vertex_succ_list
    path_dict = {i: [[i]] for i in range(len(reading_frame)) if len(reading_frame[i]) > 0}
    for i,succ_vlist in enumerate(succ_list) if gene.strand == '+' else reversed(list(enumerate(succ_list))):
        if len(succ_vlist) > 0:
            for succ_v in succ_vlist:
                if succ_v not in path_dict:
                    path_dict[succ_v] = add_v_to_list3(path_dict[i],succ_v)
                else:
                    path_dict[succ_v].extend(add_v_to_list3(path_dict[i],succ_v))
            path_dict[i] = []
            # print("--------")
            # print(i)
            # print(path_dict)
    return path_dict






