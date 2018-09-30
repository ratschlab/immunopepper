import sys
import re

# external package
import scipy as sp

# immuno module
from utils import complementary_seq,translate_dna_to_peptide

def junction_is_annotated(gene, gene_to_transcript_table, transcript_to_cds_table):
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
        #aa_str_mutated = translate_dna_to_peptide(cds_string)
        #peptide_list.append((aa_str_mutated, ts))
    gene.processed = True
    return ts_list


def find_background_peptides(gene, ref_seq, gene_to_transcript_table, transcript_cds_table):
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
        peptide_list.append((aa_str_mutated, ts))
    gene.processed = True
    return peptide_list


def peptide_match(ref_peptides, peptide):
    match_ts_list = []
    for ref in ref_peptides:
        ref_peptide = ref[0][0]
        transcript = ref[1]
        if not re.search(peptide, ref_peptide) is None:
            match_ts_list.append(transcript)
    return match_ts_list

## Returns true if there is a stop codon in the sequence. All codons that are fully in the
#  interval [start_coord<->end_coord] are checked.
# seq: Nucleotide sequence of vertex/CDS region
# start_coord: Read start coordinate
# stop_coord: Read stop coordinate
# strand: Read direction, one of {"+","-"}
def has_stop_codon_initial(seq, start_coord, end_coord, strand):
    if strand == "+":
        assert (start_coord <= end_coord)
        substr = seq[start_coord - 1:end_coord]
    else:  # strand=="-"
        assert (start_coord >= end_coord)
        substr = complementary_seq(seq[end_coord - 1:start_coord][::-1])
    for idx in range(0, len(substr) - 2, 3):
        nuc_frame = substr[idx:idx + 3]
        if nuc_frame.lower() in ["tag", "taa", "tga"]:
            return True
    return False

## Returns true if there is a stop codon found spanning the two sequences
# seq_prop: Propagating sequence
# seq_accept: Accepting sequence
# read_frame: Read frame of propagating vertex
# strand: Direction of read strand
def has_stop_codon_cross(seq_prop, seq_accept, read_frame, strand):
    if strand == "+":
        check_seq = seq_prop[-read_frame:] + seq_accept
    else:  # strand=="-"
        check_seq = seq_prop[::-1][-read_frame:] + seq_accept[::-1]
    for idx in range(0, len(check_seq) - 2, 3):
        nuc_frame = check_seq[idx:idx + 3]
        if nuc_frame.lower() in ["tag", "taa", "tga"]:
            return True
    return False


def get_true_variant_comb(variant_comb, exon_list):
    new_variant = []
    if variant_comb == ['.']:
        return '.'
    for ivariant in variant_comb:
        if int(ivariant) in range(int(exon_list[0]),int(exon_list[1])) or \
                int(ivariant) in range(int(exon_list[2]),int(exon_list[3])):
            new_variant.append(ivariant)
    return ';'.join(new_variant)


def get_exon_dict(metadata_list):
    exon_dict = {}
    for line in metadata_list:
        items = line.strip().split('\t')
        ## TODO: need to come up with a new way to index the exon list
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


def get_filtered_output_list(metadata_list,peptide_list):
    exon_list = get_exon_dict(metadata_list)
    remove_id_list = get_remove_id(exon_list)
    filtered_meta_list = []
    filtered_peptide_list = []
    for i,imeta_line in enumerate(metadata_list):
        items = imeta_line.strip().split('\t')
        idx = items[0]
        if idx not in remove_id_list:
            filtered_meta_list.append(imeta_line)
            filtered_peptide_list.append(peptide_list[i])
    return filtered_meta_list, filtered_peptide_list



