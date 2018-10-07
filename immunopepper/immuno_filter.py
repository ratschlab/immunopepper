"""Contain the function used for filtering"""
import sys
import re

import scipy as sp

from utils import complementary_seq,translate_dna_to_peptide


def junction_is_annotated(gene, gene_to_transcript_table, transcript_to_cds_table):
    """ Indicate whether exon pair also appear in transcript given by .gtf file

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
            if '%i:%i' % (vertices[1, x[0]] - 1, vertices[0, x[1]]) in transcript_junctions:
                junction_flag[x[0], x[1]] = 1
                junction_flag[x[1], x[0]] = 1

    return junction_flag


def find_background_transcript(gene, ref_seq, gene_to_transcript_table, transcript_cds_table):
    """Calculate the transcript derived from the annotation file.
    It is used for generating artificial reads.

    Parameters
    ----------
    gene: Object. Created by SplAdder
    ref_seq: List(str). Reference sequence of certain chromosome.
    gene_to_transcript_table: Dict. "Gene" -> "transcript"
    transcript_to_cds_table: Dict. "transcript" -> "cds"

    Returns
    -------
    ts_list: List[str]. List of all the peptide translated from the given
        splicegraph and annotation.
    """
    gene_transcripts = gene_to_transcript_table[gene.name]
    ts_list = []

    # Generate a background peptide for every variant transcript
    for ts in gene_transcripts:

        # No CDS entries for transcript in annotation file...
        if ts not in transcript_cds_table:
            # print("WARNING: Transcript not in CDS table")
            continue

        cds_list = transcript_cds_table[ts]

        # Reverse CDS in reverse strand transcription...
        # add 3 bases to the last cds part to account for non-annotated stops
        if gene.strand.strip() == "-":
            cds_list = cds_list[::-1]
            cds_list[-1] = (cds_list[-1][0] - 3, cds_list[-1][1], cds_list[-1][2])
        else:
            cds_list[-1] = (cds_list[-1][0], cds_list[-1][1] + 3, cds_list[-1][2])

        cds_string = ""
        first_cds = True

        # Append transcribed CDS regions to the output
        for coord_left, coord_right, frameshift in cds_list:

            # Apply initial frameshift on the first CDS of the transcript
            if first_cds:
                if gene.strand.strip() == "+":
                    coord_left += frameshift
                else:
                    coord_right -= frameshift
                first_cds = False

            nuc_seq = ref_seq[coord_left:coord_right + 1]

            # Accumulate new DNA sequence...
            if gene.strand.strip() == "+":
                cds_string += nuc_seq
            elif gene.strand.strip() == "-":
                cds_string += complementary_seq(nuc_seq[::-1])
            else:
                print("ERROR: Invalid strand...")
                sys.exit(1)
        ts_list.append(cds_string)
    gene.processed = True
    return ts_list


def find_background_peptides(gene, ref_seq, gene_to_transcript_table, transcript_cds_table):
    """Calculate the peptide translated from the complete transcript instead of exon pairs

    Parameters
    ----------
    gene: Object. Created by SplAdder
    ref_seq: List(str). Reference sequence of certain chromosome.
    gene_to_transcript_table: Dict. "Gene" -> "transcript"
    transcript_to_cds_table: Dict. "transcript" -> "cds"

    Returns
    -------
    peptide_list: List[str]. List of all the peptide translated from the given
        splicegraph and annotation.
    """
    gene_transcripts = gene_to_transcript_table[gene.name]
    peptide_list = []
    ref_seq_str = ''.join(ref_seq)
    # Generate a background peptide for every variant transcript
    for ts in gene_transcripts:

        # No CDS entries for transcript in annotation file...
        if ts not in transcript_cds_table:
            #print("WARNING: Transcript not in CDS table")
            continue

        cds_list = transcript_cds_table[ts]

        # Reverse CDS in reverse strand transcription...
        # add 3 bases to the last cds part to account for non-annotated stops
        if gene.strand.strip() == "-":
            cds_list = cds_list[::-1]
            cds_list[-1] = (cds_list[-1][0] - 3, cds_list[-1][1], cds_list[-1][2])
        else:
            cds_list[-1] = (cds_list[-1][0], cds_list[-1][1] + 3, cds_list[-1][2])

        cds_string = ""
        first_cds = True

        # Append transcribed CDS regions to the output
        for coord_left, coord_right, frameshift in cds_list:

            # Apply initial frameshift on the first CDS of the transcript
            if first_cds:
                if gene.strand.strip() == "+":
                    coord_left += frameshift
                else:
                    coord_right -= frameshift
                first_cds = False

            nuc_seq = ref_seq_str[coord_left:coord_right+1]

            # Accumulate new DNA sequence...
            if gene.strand.strip() == "+":
                cds_string += nuc_seq
            elif gene.strand.strip() == "-":
                cds_string += complementary_seq(nuc_seq[::-1])
            else:
                print("ERROR: Invalid strand...")
                sys.exit(1)

        aa_str_mutated = translate_dna_to_peptide(cds_string)
        # aa_str_mutated = (peptide, stop_codon_indicator)
        peptide_list.append(aa_str_mutated[0])
    gene.processed = True
    return peptide_list


def peptide_match(background_peptide_list, peptide):
    """ Find if the translated exon-pair peptide also appear in the annotation file.
    Parameters
    ----------
    background_peptide_list: List(str) All the peptide translated from transcripts in annotation file
    peptide: str. peptide translated from certain exon-pairs

    Returns
    -------
    count: int. Count how many matches exist between background peptide and given peptide

    """
    count = 0
    for ref in background_peptide_list:
        ref_peptide = ref
        if not re.search(peptide, ref_peptide) is None:
            count += 1
    return count


def get_true_variant_comb(variant_comb, exon_list):
    """ Filter out the variants that do not take effect

    Parameters
    ----------
    variant_comb: List(str). List of all the variant positions
    exon_list: List(str). List of the four exon terminals.

    Returns
    -------
    new_variant: List(str). List of the variants that take effect.

    """
    new_variant = []
    if variant_comb == ['.']:
        return '.'
    for ivariant in variant_comb:
        if int(ivariant) in range(int(exon_list[0]),int(exon_list[1])) or \
                int(ivariant) in range(int(exon_list[2]),int(exon_list[3])):
            new_variant.append(ivariant)
    return ';'.join(new_variant)


def get_exon_dict(metadata_list):
    """ Get an auxiliary dictionary used for filtering

    Parameters
    ----------
    metadata_list: List(str). Returned by `calculate_output_peptide`

    Returns
    -------
    exon_dict: dict. (read_frame, pos_mid_1, pos_mid_2, variant_combination)
        -> List[Tuple(output_idx, pos_start, pos_end,)]

    """

    exon_dict = {}
    for line in metadata_list:
        items = line.strip().split('\t')
        #(Jiayu) TODO: need to come up with a new way to index the exon list
        exon_list = items[-4].split(';')
        idx = items[0]
        read_frame = items[1]
        strand = items[4]
        variant_comb = items[-6].split(';')
        variant_comb = get_true_variant_comb(variant_comb,exon_list)
        if strand == '+':
            key = (read_frame,exon_list[1],exon_list[2],variant_comb)
            if key in exon_dict:
                exon_dict[key].append((idx, exon_list[0],exon_list[3]))
            else:
                exon_dict[key] = [(idx, exon_list[0],exon_list[3])]
        else:
            key = (read_frame,exon_list[3], exon_list[0],variant_comb)
            if key in exon_dict:
                exon_dict[key].append((idx, exon_list[2], exon_list[1]))
            else:
                exon_dict[key] = [(idx, exon_list[2], exon_list[1])]
    return exon_dict


def get_remove_id(metadata_dict):
    """ Get all redundant output id.
    Parameters
    ----------
    metadata_dict: Dict. Returned by get_exon_dict

    Returns
    -------
    remove_id_list: List[str]. The list of id to be removed.

    """
    remove_id_list = []
    for exon_pair_list in metadata_dict.values():
        L = len(exon_pair_list)
        if L < 2:
            continue
        for i, exon_pair in enumerate(exon_pair_list):
            for j in range(L):
                output_id = exon_pair[0]
                i_pos1 = exon_pair[1]
                i_pos2 = exon_pair[2]
                j_pos1 = exon_pair_list[j][1]
                j_pos2 = exon_pair_list[j][2]
                if j!=i and i_pos1 >= j_pos1 and i_pos2 <= j_pos2:
                    remove_id_list.append(output_id)
                    break
    return remove_id_list


def get_filtered_output_list(metadata_list, peptide_list, expr_lists):
    """ Filter redundant output

    Parameters
    ----------
    output_peptide_list: List[str]. Contain all the possible output peptide in the given splicegraph.
    output_metadata_list: List[str]. Contain the correpsonding medata data for each output peptide.
    expr_list: List[List[Tuple]]. The Tuple is (len_segments, expression). The Inner List is all the component
        segment corresponding to one exon pair. The outer List is the aggregation for all the exon-pair segment expression
        information in the given gene.
    All the above three are returned from `calculate_output_peptide`
    """
    exon_list = get_exon_dict(metadata_list)
    remove_id_list = get_remove_id(exon_list)
    filtered_meta_list = []
    filtered_peptide_list = []
    filtered_expr_lists = []
    for i,imeta_line in enumerate(metadata_list):
        items = imeta_line.strip().split('\t')
        idx = items[0]
        if idx not in remove_id_list:
            filtered_meta_list.append(imeta_line)
            filtered_peptide_list.append(peptide_list[i])
            filtered_expr_lists.append(expr_lists[i])
    return filtered_meta_list, filtered_peptide_list, filtered_expr_lists



