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
                
            nuc_seq = ref_seq[coord_left:coord_right+1]

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
