"""Countains code related to translation"""

from .utils import get_exon_expr,get_sub_mut_dna
from .immuno_nametuple import Coord
from .immuno_nametuple import Flag
from .immuno_nametuple import Peptide
from .immuno_nametuple import ReadingFrameTuple
from .constant import NOT_EXIST

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


def get_full_peptide(gene, seq, cds_list, countinfo, Idx, mode):
    """ 
    Output translated peptide and segment expression list given cds_list

    Parameters
    ----------
    gene: Object, created by SplAdder.
    seq: str. Gene sequence.
    cds_list: List[Tuple(v_start,v_stop,reading_frame)]
    countinfo: Namedtuple, SplAdder count information
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
        cds_expr = get_exon_expr(gene, coord_left, coord_right, countinfo, Idx)
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


def complementary_seq(dna_seq):
    """ Yields the complementary DNA sequence
    Only convert the character in comp_dict.
    Otherwise remain the same.
    """
    comp_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([comp_dict[nuc] if nuc in comp_dict else nuc for nuc in dna_seq])


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
        peptide_dna_str_mut = get_sub_mut_dna(mut_seq, modi_coord, variant_comb, somatic_mutation_sub_dict, strand)
        peptide_dna_str_ref = get_sub_mut_dna(ref_seq, modi_coord, NOT_EXIST, somatic_mutation_sub_dict, strand)
    else:  # strand == "-"
        peptide_dna_str_mut = complementary_seq(get_sub_mut_dna(mut_seq, modi_coord, variant_comb, somatic_mutation_sub_dict, strand))
        peptide_dna_str_ref = complementary_seq(get_sub_mut_dna(ref_seq, modi_coord, NOT_EXIST, somatic_mutation_sub_dict, strand))
    peptide_mut, mut_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_mut)
    peptide_ref, ref_has_stop_codon = translate_dna_to_peptide(peptide_dna_str_ref)

    # if the stop codon appears before translating the second exon, mark 'single'
    if modi_coord.start_v2 == NOT_EXIST or len(peptide_mut)*3 <= abs(modi_coord.stop_v1 - modi_coord.start_v1) + 1:
        is_isolated = True
    else:
        is_isolated = False
    peptide = Peptide(peptide_mut, peptide_ref)
    flag = Flag(mut_has_stop_codon, is_isolated)
    return peptide,flag


