"""Contain the function used for filtering"""
import sys
import re

import scipy as sp

from utils import complementary_seq,translate_dna_to_peptide,get_exon_expr


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
            print("ERROR: Invalid strand...")
            sys.exit(1)
    cds_peptide, is_stop_flag = translate_dna_to_peptide(cds_string)
    return cds_expr_list, cds_string, cds_peptide


def get_full_peptide_list(gene, ref_seq, all_path_dict, Segments, Idx):
    """
    Output all peptides and corresponding expression list given all possible paths
    """
    def re_arrange_dict(all_path_dict):
        """Index path with start vertices"""
        new_dict = {}
        for path_list in all_path_dict.values():
            for path in path_list:
                if path[0] not in new_dict:
                    new_dict[path[0]] = [path]
                else:
                    new_dict[path[0]].append(path)
        return new_dict

    def convert_to_cds_list(new_dict,gene):
        """Change vertex index to cds tuple"""
        all_cds_list = []
        vlist = gene.splicegraph.vertices
        for start_point,path_list in new_dict.items():
            reading_frame_list = list(gene.splicegraph.reading_frames[start_point])
            for reading_frame_tuple in reading_frame_list:
                for path in path_list:
                    cds_list = []
                    cds_list.append(reading_frame_tuple)
                    for succ_point in path[1:]:
                        cds_list.append((vlist[0][succ_point],vlist[1][succ_point],0))
                    all_cds_list.append(cds_list)
        return all_cds_list

    def convert_cds_list_to_string(cds_list):
        """For the header line in the output"""
        cds_str_list = [str(cds_item) for cds_item in cds_list]
        return '_'.join(cds_str_list)

    new_path_dict = re_arrange_dict(all_path_dict)
    all_cds_list = convert_to_cds_list(new_path_dict,gene)
    peptide_list = []
    expr_lists = []
    for i,cds_list in enumerate(all_cds_list):
        cds_expr_list, cds_string, cds_peptide = get_full_peptide(gene, ref_seq, cds_list, Segments, Idx, mode='full')
        peptide_list.append((convert_cds_list_to_string(cds_list),cds_peptide))
        expr_lists.append(cds_expr_list)
    return peptide_list, expr_lists


def find_background_peptides(gene, ref_seq, gene_to_transcript_table, transcript_cds_table, Segments, Idx):
    """Calculate the peptide translated from the complete transcript instead of single exon pairs

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
    (ts_list): List[str]. List of all the transcript indicated by the  annotation file
        can be used to generate artifical reads.
    """
    gene_transcripts = gene_to_transcript_table[gene.name]
    peptide_list = []
    expr_lists = []
    # Generate a background peptide for every variant transcript
    for ts in gene_transcripts:

        # No CDS entries for transcript in annotation file...
        if ts not in transcript_cds_table:
            #print("WARNING: Transcript not in CDS table")
            continue
        cds_list = transcript_cds_table[ts]
        cds_expr_list, cds_string, cds_peptide = get_full_peptide(gene,ref_seq,cds_list,Segments,Idx,mode='back')
        peptide_list.append((ts,cds_peptide))
        expr_lists.append(cds_expr_list)
    gene.processed = True
    return peptide_list, expr_lists


def peptide_match(ref_peptides, peptide):
    """ Find if the translated exon-pair peptide also appear in the annotation file.

    Parameters
    ----------
    background_peptide_list: List(str) All the peptide translated from transcripts in annotation file
    peptide: str. peptide translated from certain exon-pairs

    Returns
    -------
    count: int. Count how many matches exist between background peptide and given peptide
    """
    match_ts_list = []
    for ref in ref_peptides:
        ref_peptide = ref[1]
        transcript = ref[0]
        if not re.search(peptide, ref_peptide) is None:
            match_ts_list.append(transcript)
    return match_ts_list


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
        ## TODO: need to come up with a new way to index the exon list
        exon_list = items[-4].split(';')
        idx = items[0]
        read_frame = items[1]
        strand = items[4]
        if strand == '+':
            key = (read_frame,exon_list[1],exon_list[2])
            if key in exon_dict:
                exon_dict[key].append((idx, exon_list[0],exon_list[3]))
            else:
                exon_dict[key] = [(idx, exon_list[0],exon_list[3])]
        else:
            key = (read_frame,exon_list[3], exon_list[0])
            if key in exon_dict:
                exon_dict[key].append((idx, exon_list[2], exon_list[1]))
            else:
                exon_dict[key] = [(idx, exon_list[2], exon_list[1])]
    return exon_dict


def get_remove_id(metadata_dict):
    # TODO: [Warning!] if two output lines have identical
    # coordinates and readframe
    # both of them will be removed. Need to fix.
    remove_id_list = []
    for exon_pair_list in metadata_dict.values():
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


def get_filtered_output_list(metadata_list,peptide_list,expr_lists):
    """ Get all redundant output id.

    Parameters
    ----------
    metadata_dict: Dict. Returned by get_exon_dict

    Returns
    -------
    remove_id_list: List[str]. The list of id to be removed.
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
    return filtered_meta_list, filtered_peptide_list,filtered_expr_lists



