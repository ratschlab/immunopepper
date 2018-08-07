import itertools

import scipy as sp
import numpy as np

# Converts a binary adjacency matrix to a list of directed edges
def to_adj_list(adj_matrix):
    adj_list = []
    assert (adj_matrix.shape[0] == adj_matrix.shape[1])
    for idx in range(adj_matrix.shape[0]):
        for jdx in range(adj_matrix.shape[0]):
            if adj_matrix[idx, jdx] == 1 and idx <= jdx:
                adj_list.append([idx, jdx])
    return adj_list

# Returns a list of successors by vertex, sensitive to the read strand
def to_adj_succ_list(adj_matrix, vertex_map, read_strand):
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


# Find overlapping CDS given a list of CDS starts
def find_overlapping_cds_simple(v_start, v_stop, cds_begins, strand):
    if strand == "+":
        return filter(lambda cds_begin: cds_begin[0] >= v_start and cds_begin[1] <= v_stop, cds_begins)
    else:
        return filter(lambda cds_begin: cds_begin[1] >= v_start and cds_begin[0] <= v_stop, cds_begins)


def attribute_list_to_dict(a_list):
    a_dict = {}
    for attribute_pair in a_list:
        pair = attribute_pair.split(' ')
        a_dict[pair[0]] = pair[1][1:-1]  # delete "", currently now work on level 2
    return a_dict


# Returns strand-sensitive order between two genomic coordinates
def leq_strand(coord1, coord2, strand_mode):
    if strand_mode == "+":
        return coord1 <= coord2
    else:  # strand_mode == "-"
        return coord1 >= coord2

# Yields the complementary DNA sequence
# dna_seq: Input nucleotide sequence
def complementary_seq(dna_seq):
    comp_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    comp_dict_keys = comp_dict.keys()
    return "".join(map(lambda nuc: comp_dict[nuc] if nuc in comp_dict_keys else nuc, dna_seq))


# Returns header labels corresponding to donor_id and mutation_mode
def header_labels(donor_id, mutation_mode):
    if mutation_mode is None:
        mutation_type = "REFERENCE"
    elif mutation_mode == "both":
        mutation_type = "GERM_SOMATIC"
    elif mutation_mode == "somatic_only":
        mutation_type = "SOMATIC"
    elif mutation_mode == "germline_only":
        mutation_type = "GERM"

    peptide_type = "REFERENCE" if donor_id == "ref" else donor_id

    return (peptide_type, mutation_type)

# Encodes chromosome to same cn
def encode_chromosome(in_num):
    convert_dict = {23: "X", 24: "Y", 25: "MT"}
    return convert_dict[in_num] if in_num in convert_dict.keys() else str(in_num)

# Translate a DNA sequence encoding a peptide to amino-acid sequence via RNA
def translate_dna_to_peptide(dna_str):
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
    if start_v2 != '.':
        if strand == '-':
            sub_dna_list = list(background_seq[start_v1:stop_v1][::-1] + background_seq[start_v2:stop_v2][::-1])
        else:
            sub_dna_list = list(background_seq[start_v1:stop_v1] + background_seq[start_v2:stop_v2])
    else:
        if strand == '-':
            sub_dna_list = list(background_seq[start_v1:stop_v1][::-1])
        else:
            sub_dna_list = list(background_seq[start_v1:stop_v1])
    if len(variant_comb) == 1: #no mutation exist
        return ''.join(sub_dna_list)
    for variant_ipos in variant_comb:
        ref_base = mutation_sub_dic_maf[variant_ipos]['ref_base']
        mut_base = mutation_sub_dic_maf[variant_ipos]['mut_base']
        strand = mutation_sub_dic_maf[variant_ipos]['strand']
        variant_Classification = mutation_sub_dic_maf[variant_ipos]['variant_Classification']
        variant_Type = mutation_sub_dic_maf[variant_ipos]['variant_Type']
        sub_dna_list[variant_ipos-start_v1] = mut_base
    return ''.join(sub_dna_list) 


def cross_peptide_result(read_frame, strand, variant_comb, mutation_sub_dic_maf, background_seq, peptide_accept_coord):
    """

    Parameters
    ----------
    read_frame: tuple, (read_start_codon, read_stop_codon, emitting_frame)
    strand: str, '+' or '-'
    mut_seq: str, mutation string
    background_seq: str, background string
    peptide_accept_coord: np.ndarray, the 'next' vertex coordinate

    Returns
    -------
    peptide_mut: str, the translated peptide of mutation sequence
    peptide_ref: str, the translated peptide of reference sequence
    start_v1, stop_v1: int, the corresponding start and end position of sequence of first exon which can output peptide
    start_v2, stop_v2:
    mut_has_stop_codon: bool, flag indicating if the sequence of interest has stop codon
    is_single: bool, flag indicating if the outputed peptide only comes from one exon
    next_reading_frame: tuple, the reading frame to be propagated to the next vertex

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


    # python is 0-based while gene annotation file(.gtf) is one based
    # so we need to do a little modification
    if strand == "+":
        start_v2 = peptide_accept_coord[0]
        stop_v2 = peptide_accept_coord[1] - next_emitting_frame
        peptide_dna_str_mut = get_sub_mut_dna(background_seq, start_v1, stop_v1, start_v2, stop_v2, variant_comb, mutation_sub_dic_maf, strand)
        peptide_dna_str_ref = background_seq[start_v1:stop_v1] + background_seq[start_v2:stop_v2]
        next_start_v1 = start_v2 + accepting_frame
        next_stop_v1 = peptide_accept_coord[1]
    else:  # strand == "-"
        start_v2 = peptide_accept_coord[0] + next_emitting_frame
        stop_v2 = peptide_accept_coord[1]
        peptide_dna_str_mut = complementary_seq(get_sub_mut_dna(background_seq, start_v1, stop_v1, start_v2, stop_v2, variant_comb, mutation_sub_dic_maf, strand))
        peptide_dna_str_ref = complementary_seq(background_seq[start_v1:stop_v1][::-1] + background_seq[start_v2:stop_v2][::-1])
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

    return peptide_mut, peptide_ref, start_v1, stop_v1, start_v2, stop_v2, \
           mut_has_stop_codon,is_isolated,next_reading_frame


def isolated_peptide_result(read_frame, strand, variant_comb, mutation_sub_dic_maf,background_seq):
    """
    Deal with translating isolated cds, almost the same as cross_peptide_result

    Parameters
    ----------
    read_frame
    strand
    mut_seq
    background_seq

    Returns
    -------
    is_isolated: Bool, True as Default
    mut_has_stop_codon: Bool, assert to be True
    start_v2, stop_v2: '-' means not exist

    """
    start_v1, stop_v1, emitting_frame = read_frame
    start_v2 = '.'  # does not exist
    stop_v2 = '.'  # does not exist
    if strand == '+':
        peptide_dna_str_mut = get_sub_mut_dna(background_seq, start_v1, stop_v1, start_v2, stop_v2, variant_comb, mutation_sub_dic_maf, strand)
        peptide_dna_str_ref = background_seq[start_v1:stop_v1]
    else:
        peptide_dna_str_mut = complementary_seq(get_sub_mut_dna(background_seq, start_v1, stop_v1, start_v2, stop_v2, variant_comb, mutation_sub_dic_maf, strand))
        peptide_dna_str_ref = background_seq[start_v1:stop_v1][::-1]

    peptide_mut, mut_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_mut)
    peptide_ref, ref_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_ref)

    is_isolated = True

    return peptide_mut, peptide_ref, start_v1, stop_v1, start_v2, stop_v2, mut_has_stop_codon, is_isolated


def get_size_factor(strains,lib_file_path):
    libs = sp.loadtxt(lib_file_path, dtype='str', skiprows=1, delimiter='\t')
    a, b = sp.where(strains[:, sp.newaxis] == libs[:, 0])
    assert sp.all(libs[b, 0] == strains)
    libs = libs[b, :]
    med = sp.median(libs[:, 1].astype('float'))
    sf = med / libs[:, 1].astype('float')
    return sf

def get_all_comb(array,r=None):
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
    If a vertex has a successor, it is not isolated.

    Generally if it has read frame, this read frame should be propogated by other
    vertex. Therefore this vertex should not be isolated, except the very beginning
    vertex, where the reading frame is assigned by gtf file. We mark all vertices that have
    no incoming or outgoing edges as isolated.

    Parameters
    ----------
    gene: gene object
    idx: the id of vertex

    Returns
    -------

    """

    if len(gene.vertex_succ_list[idx]) > 0:
        return False

    return sp.sum(gene.splicegraph.edges[:, idx]) == 0
    #sg = gene.splicegraph
    ## AK TODO: this needs fixing. pop() will remove the reading frame, also it will pick an arbitrary one.
    #coord_diff = sg.vertices[0, idx] - sg.reading_frames[idx].pop()[0]
    #print(coord_diff)
    #if abs(coord_diff) > 10:
    #    return True

def is_output_redundant(gene, start_v1, stop_v1, start_v2, stop_v2):
    strand = gene.strand
    output_vertex_dict = gene.output_vertex_dict
    if strand == '+':
        if (stop_v1,start_v2) in output_vertex_dict:
            all_output_comb = output_vertex_dict[(stop_v1,start_v2)]
            for comb in all_output_comb:
                if start_v1 >= comb[0] and stop_v2 <= comb[1]:
                    return True
            output_vertex_dict[(stop_v1,start_v2)].append((start_v1,stop_v2))
        else:
            output_vertex_dict[(stop_v1,start_v2)] = [(start_v1,stop_v2)]
    else:  # strand == '-'
        if (start_v1, stop_v2) in output_vertex_dict:
            all_output_comb = output_vertex_dict[(start_v1, stop_v2)]
            for comb in all_output_comb:
                if stop_v1 <= comb[0] and start_v2 >= comb[1]:
                    return True
            output_vertex_dict[(start_v1, stop_v2)].append((stop_v1, start_v2))
        else:
            output_vertex_dict[(start_v1, stop_v2)] = [(stop_v1, start_v2)]
    return False

def is_in_junction_list(v1,v2,strand,junction_list):

    return int(':'.join([str(v1[1]),  str(v2[0]), strand]) in junction_list)

   # coordi_match = np.logical_and(junction_list[:,0] == str(v1[1]),junction_list[:,1] == str(v2[0]))
   # final_match  = sum(np.logical_and(coordi_match,junction_list[:,2]==strand))
   # return final_match
