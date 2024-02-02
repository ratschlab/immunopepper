"""Contains all the output computation based on gene splicegraph"""
from collections import defaultdict
import numpy as np
import logging
import timeit

from immunopepper.filter import filter_redundant_junctions
from immunopepper.filter import junctions_annotated
from immunopepper.filter import is_intron_in_junction_list
from immunopepper.io_ import namedtuple_to_str
from immunopepper.io_ import save_fg_peptide_set
from immunopepper.io_ import save_kmer_matrix
from immunopepper.mutations import get_mutated_sequence
from immunopepper.mutations import exon_to_mutations
from immunopepper.mutations import get_mut_comb
from immunopepper.mutations import exon_to_expression
from immunopepper.namedtuples import Coord
from immunopepper.namedtuples import OutputPeptide
from immunopepper.namedtuples import OutputMetadata
from immunopepper.namedtuples import ReadingFrameTuple
from immunopepper.namedtuples import VertexPair
from immunopepper.preprocess import search_edge_metadata_segmentgraph
from immunopepper.translate import get_exhaustive_reading_frames
from immunopepper.translate import get_full_peptide
from immunopepper.translate import isolated_peptide_result
from immunopepper.translate import get_peptide_result
from immunopepper.translate import cross_peptide_result
from immunopepper.utils import get_segment_expr
from immunopepper.utils import replace_I_with_L


def collect_background_transcripts(gene=None, ref_seq_file=None, chrm=None, mutation=None):
    """Calculate the background peptide""

       Parameters
       ----------
       gene: Object, returned by SplAdder.
       ref_seq: Str, reference sequnce of specific chromosome
       mutation: Namedtuple Mutation, store the mutation information of specific chromosome and sample.
           has the attribute ['mode', 'maf_dict', 'vcf_dict']

       Returns
       -------
       ref_mut_seq: Dict. (sequence_type) -> list[char].
       """
    gene.from_sparse()
    min_pos = gene.splicegraph.vertices.min()
    max_pos = gene.splicegraph.vertices.max()

    # apply germline mutation
    # when germline mutation is applied, background_seq != ref_seq
    # otherwise, background_seq = ref_seq
    ref_mut_seq = get_mutated_sequence(fasta_file=ref_seq_file,
                                          chromosome=chrm,
                                          pos_start=min_pos,
                                          pos_end=max_pos,
                                          mutation_dict=mutation.germline_dict)
    return ref_mut_seq



def collect_vertex_pairs(gene=None, gene_info=None, ref_seq_file=None, chrm=None, idx=None, mutation=None, all_read_frames=False, disable_concat=False, kmer_length=None, filter_redundant=False):
    """Calculate the output peptide for every exon-pair in the splicegraph

       Parameters
       ----------
       gene: Object, returned by SplAdder.
       ref_seq: Str, reference sequnce of specific chromosome
       idx: Namedtuple Idx, has attribute idx.gene and idx.sample
       mutation: Namedtuple Mutation, store the mutation information of specific chromosome and sample.
           has the attribute ['mode', 'maf_dict', 'vcf_dict']
       disable_concat: bool, flag indicating whether to disable the concanation of vertices into triples
       kmer_length: int. kmer length
       filter_redundant: flag indicating whether to remove pairs spanning the same intron

       Returns
       -------
       vertex_pairs: List of VertexPair.
       ref_mut_seq: Dict. (sequence_type) -> list[char].
       exon_som_dict: Dict. (exon_id) |-> (mutation_postion)
       """

    gene.from_sparse()
    sg = gene.splicegraph
    min_pos = gene.splicegraph.vertices.min()
    max_pos = gene.splicegraph.vertices.max()

    output_id = 0

    # apply germline mutation
    # when germline mutation is applied, background_seq != ref_seq
    # otherwise, background_seq = ref_seq
    ref_mut_seq = get_mutated_sequence(fasta_file=ref_seq_file,
                                          chromosome=chrm,
                                          pos_start=min_pos,
                                          pos_end=max_pos,
                                          mutation_dict=mutation.germline_dict)

    # apply somatic mutation
    # exon_som_dict: (exon_id) |-> (mutation_postion)
    exon_som_dict = None
    if mutation.somatic_dict is not None:
        exon_som_dict = exon_to_mutations(gene, mutation.somatic_dict)

    vertex_pair_list = []
    if all_read_frames:
        reading_frame_dict = dict(get_exhaustive_reading_frames(sg, gene.strand, gene_info.vertex_order))
    else: # use reading frames from annotation
        reading_frame_dict = dict(gene_info.reading_frames)

    for v_id in gene_info.vertex_order:
        n_read_frames = len(reading_frame_dict[v_id])
        if n_read_frames == 0:  # no cds start, skip the vertex
            continue
        if len(gene_info.vertex_succ_list[v_id]) == 0: # if no successive vertex, we set it to np.nan, translate and output it
            gene_info.vertex_succ_list[v_id].append(np.nan)

        for prop_vertex in gene_info.vertex_succ_list[v_id]:
            vertex_list = [v_id, prop_vertex]
            mut_seq_comb = get_mut_comb(exon_som_dict,vertex_list)
            for read_frame_tuple in sorted(reading_frame_dict[v_id]):
                has_stop_flag = True
                for variant_comb in mut_seq_comb:  # go through each variant combination
                    if prop_vertex is not np.nan:
                        peptide, modi_coord, flag, next_start_v1, next_stop_v1, next_emitting_frame = cross_peptide_result(read_frame_tuple, gene.strand, variant_comb, mutation.somatic_dict, ref_mut_seq, sg.vertices[:, prop_vertex], min_pos, all_read_frames)
                        orig_coord = Coord(sg.vertices[0, v_id], sg.vertices[1, v_id], sg.vertices[0, prop_vertex], sg.vertices[1, prop_vertex])
                        # no propagation needed in all reading frame mode,  RF not labelled as annotated
                        if (not flag.has_stop) and (not all_read_frames) \
                            and (not ReadingFrameTuple(next_start_v1, next_stop_v1,
                                                       next_emitting_frame, True) in reading_frame_dict[prop_vertex] ):
                            reading_frame_dict[prop_vertex].add(ReadingFrameTuple(next_start_v1,
                                                                                  next_stop_v1, next_emitting_frame, False))
                    else:
                        peptide, modi_coord, flag = isolated_peptide_result(read_frame_tuple, gene.strand, variant_comb, mutation.somatic_dict, ref_mut_seq, min_pos, all_read_frames)
                        orig_coord = Coord(sg.vertices[0, v_id],sg.vertices[1, v_id], np.nan, np.nan)
                    has_stop_flag = has_stop_flag and flag.has_stop
                gene_outputid = str(idx.gene) + ':' + str(output_id)
                vertex_pair = VertexPair(output_id=gene_outputid,
                                         read_frame=read_frame_tuple,
                                         modified_exons_coord=modi_coord,
                                         original_exons_coord=orig_coord,
                                         vertex_idxs=vertex_list,
                                         has_stop_codon=has_stop_flag,
                                         peptide_weight="{:.3f}".format(1/n_read_frames))
                vertex_pair_list.append(vertex_pair)
                output_id += 1

    if filter_redundant:
        vertex_pair_list = filter_redundant_junctions(vertex_pair_list, gene.strand)
    concat_vertex_pair_list = defaultdict(list, {'2-exons': vertex_pair_list})
    if not disable_concat:
        concat_vertex_pair_list['3-exons'] = collect_vertex_triples(gene, vertex_pair_list, kmer_length)
    #vertex_pair_list += concat_vertex_pair_list

    return concat_vertex_pair_list, ref_mut_seq, exon_som_dict


def collect_vertex_triples(gene, vertex_pairs, k):
    """

    Parameters
    ----------
    gene: Object, returned by SplAdder.
    vertex_pairs: List of VertexPair.
    k: Int. the length of kmers.

    Returns
    -------
    concat_simple_meta_list: List of VertexPair, specifically for triple-vertice cases.

    """

    def _in_the_same_read_frame(front_coord, back_coord, strand):
        if strand == '+':
            return (front_coord.stop_v2 - back_coord.start_v1) % 3 == 0
        else:
            return (back_coord.stop_v1 - front_coord.start_v2) % 3 == 0

    concat_vertex_pair_list = []
    vertex_lens = gene.splicegraph.vertices[1,:] - gene.splicegraph.vertices[0,:]
    key_id_list = np.where(vertex_lens < (k + 1) * 3)[0]
    for key_id in key_id_list:
        front_id_list = [i for i, vp in enumerate(vertex_pairs) if vp.vertex_idxs[1] == key_id
                         and not vp.has_stop_codon]
        back_id_list = [i for i, vp in enumerate(vertex_pairs) if vp.vertex_idxs[0] == key_id
                        and vp.vertex_idxs[1] is not np.nan]
        for front_id in front_id_list:
            for back_id in back_id_list:
                front_pair = vertex_pairs[front_id]
                back_pair = vertex_pairs[back_id]
                if not _in_the_same_read_frame(front_pair.modified_exons_coord, back_pair.modified_exons_coord, gene.strand):
                    continue
                middle_exon_coord = gene.splicegraph.vertices[:, front_pair.vertex_idxs[1]]
                triple_modi_coord = Coord(start_v1=front_pair.modified_exons_coord.start_v1,
                                          stop_v1=front_pair.modified_exons_coord.stop_v1,
                                          start_v2=middle_exon_coord[0],
                                          stop_v2=middle_exon_coord[1],
                                          start_v3=back_pair.modified_exons_coord.start_v2,
                                          stop_v3=back_pair.modified_exons_coord.stop_v2)
                triple_orig_coord = Coord(start_v1=front_pair.original_exons_coord.start_v1,
                                          stop_v1=front_pair.original_exons_coord.stop_v1,
                                          start_v2=middle_exon_coord[0],
                                          stop_v2=middle_exon_coord[1],
                                          start_v3=back_pair.original_exons_coord.start_v2,
                                          stop_v3=back_pair.original_exons_coord.stop_v2)
                triple_output_id = front_pair.output_id + '_' + back_pair.output_id.split('.')[-1]
                triple_vertex_idxs = front_pair.vertex_idxs + [back_pair.vertex_idxs[-1]]
                new_vertex_triple = VertexPair(output_id=triple_output_id,
                                             read_frame=front_pair.read_frame,
                                             modified_exons_coord=triple_modi_coord,
                                             original_exons_coord=triple_orig_coord,
                                             vertex_idxs=triple_vertex_idxs,
                                             has_stop_codon=back_pair.has_stop_codon,
                                             peptide_weight=front_pair.peptide_weight)
                concat_vertex_pair_list.append(new_vertex_triple)

    return concat_vertex_pair_list


def get_and_write_peptide_and_kmer(peptide_set: object = None,
                                   gene: object = None, all_vertex_pairs: object = None, ref_mut_seq: object = None, idx: object = None,
                                   exon_som_dict: object = None, countinfo: object = None,
                                   edge_idxs: object = None, edge_counts: object = None, seg_counts: object = None,
                                   mutation: object = None, mut_count_id: object = None, table: object = None,
                                   size_factor: object = None, junction_list: object = None, kmer_database: object = None,
                                   filepointer: object = None,
                                   force_ref_peptides: object = False, kmer: object = None,
                                   all_read_frames: object = None, graph_output_samples_ids: object = None,
                                   graph_samples: object = None, out_dir: object = None, verbose_save: object = None, 
                                   fasta_save: object = None) -> object:
    """

    Parameters
    ----------
    peptide_set: set(OutputMetadata, OutputMetadata) with OutputMetadata namedtuple
    gene: Object, returned by SplAdder.
    all_vertex_pairs: List of VertexPair
    ref_mut_seq: Str, reference sequnce of specific chromosome
    idx: Namedtuple Idx, has attribute idx.gene and idx.sample
    exon_som_dict: Dict. (exon_id) |-> (mutation_postion)
    countinfo: Namedtuple, contains SplAdder count information
    mutation: Namedtuple Mutation, store the mutation information of specific chromosome and sample.
        has the attribute ['mode', 'maf_dict', 'vcf_dict']
    table: Namedtuple GeneTable, store the gene-transcript-cds mapping tables derived
       from .gtf file. has attribute ['gene_to_cds_begin', 'ts_to_cds', 'gene_to_cds']
    size_factor: Scalar. To adjust the expression counts based on the external file `libsize.tsv`
    junction_list: List. Work as a filter to indicate some exon pair has certain
       ordinary intron which can be ignored further.
    kmer_database: Set. kmers to be removed on the fly from the kmer sample or matrix files
    filepointer: namedtuple, contains the columns and paths of each file of interest
    force_ref_peptides: bool, flag indicating whether to force output of
        mutated peptides which are the same as reference peptides
    kmer: list containing the length of the kmers requested
    out_dir: str, base direactory used for temporary files
    graph_samples: list, samples contained in the splicing graph object
    fasta_save: bool. whether to save a fasta file with the peptides
    """
    # check whether the junction (specific combination of vertices) also is annotated
    # as a  junction of a protein coding transcript
    len_pep_save = 9999
    gene_annot_jx = junctions_annotated(gene, table.gene_to_ts, table.ts_to_cds)
    som_exp_dict = exon_to_expression(gene, list(mutation.somatic_dict.keys()), countinfo, seg_counts, mut_count_id)
    gene_kmer_coord = set()

    ### iterate over all vertex pairs and translate
    for kmer_type, vertex_pairs in all_vertex_pairs.items():
        time_stamp = timeit.default_timer() 
        for ii,vertex_pair in enumerate(vertex_pairs):
            new_time = timeit.default_timer()
            print(f'working on vertex pair {ii} out of {len(vertex_pairs)} - last took {new_time - time_stamp}') # AK
            time_stamp = new_time
            modi_coord = vertex_pair.modified_exons_coord
            vertex_list = vertex_pair.vertex_idxs
            tran_start_pos = modi_coord.start_v1 if gene.strand == '+' else modi_coord.stop_v1
            mut_seq_comb = get_mut_comb(exon_som_dict, vertex_pair.vertex_idxs)

            variant_id = 0
            for variant_comb in mut_seq_comb:  # go through each variant combination
                peptide,flag = get_peptide_result(vertex_pair, gene.strand, variant_comb, mutation.somatic_dict, ref_mut_seq, np.min(gene.splicegraph.vertices), all_read_frames)

                # If cross junction peptide has a stop-codon in it, the frame
                # will not be propagated because the read is truncated before it reaches the end of the exon.
                # also in mutation mode, only output the case where ref is different from mutated


                for pep_idx in np.arange(len(peptide.mut)):
                    # Do not output peptide if:
                    # (1) peptide is empty peptide
                    # (2) In mutation mode the peptide is the same as the reference, unless the user forced redundancy
                    if not peptide.mut[pep_idx] \
                            or ((mutation.mode != 'ref') and (peptide.mut[pep_idx] in peptide.ref) and (not force_ref_peptides)):
                        continue

                    # collect flags
                    new_output_id = ':'.join([gene.name, '_'.join([str(v) for v in vertex_list]), str(variant_id), str(tran_start_pos), kmer_type])
                    is_intron_in_junction_list_flag = is_intron_in_junction_list(gene.splicegraph, vertex_list, gene.strand, junction_list)

                    # collect mutations
                    if variant_comb is not np.nan and som_exp_dict is not None:  # which means there exist mutations
                        seg_exp_variant_comb = [int(som_exp_dict[ipos]) for ipos in variant_comb]
                    else:
                        seg_exp_variant_comb = np.nan  # if no mutation or no count file,  the segment expression is .


                    ### Peptides
                    peptide_set.add(namedtuple_to_str(OutputMetadata(peptide=peptide.mut[pep_idx],
                                       output_id=new_output_id,
                                       read_frame=vertex_pair.read_frame.read_phase,
                                       read_frame_annotated=vertex_pair.read_frame.annotated_RF,
                                       gene_name=gene.name,
                                       gene_chr=gene.chr,
                                       gene_strand=gene.strand,
                                       mutation_mode=mutation.mode,
                                       has_stop_codon=int(flag.has_stop),
                                       is_in_junction_list=is_intron_in_junction_list_flag,
                                       is_isolated=int(flag.is_isolated),
                                       variant_comb=variant_comb,
                                       variant_seg_expr=seg_exp_variant_comb,
                                       modified_exons_coord=modi_coord,
                                       original_exons_coord=vertex_pair.original_exons_coord,
                                       vertex_idx=vertex_list,
                                       kmer_type=kmer_type
                                       ), sep = '\t'))
                    variant_id += 1
                    output_peptide = OutputPeptide(output_id=new_output_id,
                                                    peptide=peptide.mut[pep_idx],
                                                    exons_coor=modi_coord,
                                                    strand=gene.strand,
                                                    read_frame_annotated=vertex_pair.read_frame.annotated_RF)

                    ### kmers
                    create_kmer(output_peptide, kmer, gene_kmer_coord, kmer_database)

                    if len(peptide_set) > len_pep_save:
                        save_fg_peptide_set(peptide_set, filepointer, out_dir, fasta_save,
                                            verbose=False, id=f'{kmer_type}{ii}')
                        peptide_set.clear()

        if not gene.splicegraph.edges is None:
            gene.to_sparse()

    prepare_output_kmer(gene, idx, countinfo, seg_counts, edge_idxs, edge_counts,
                        gene_kmer_coord, gene_annot_jx,
                        graph_output_samples_ids,
                        graph_samples, filepointer, out_dir, verbose=verbose_save)

    save_fg_peptide_set(peptide_set, filepointer, out_dir, fasta_save,
                        verbose=False, id=f'{kmer_type}{ii}')

def get_spanning_index(coord, k):
    """
    Generate the indexes at which a given k-mer crosses the first, second or both junctions
    :param coord: coordinate namedtuple
    :param k: length of the k-mer
    :return: spanning_id_range1, spanning_id_range2, spanning_id_range1_2: set containing all
    the start positions of the k-mer inside the peptide, such that the k-mer crosses the first junction,
    the second junction or both respectively
    """
    L1 = coord.stop_v1-coord.start_v1
    if coord.start_v2 is np.nan:
        return [np.nan], [np.nan], [np.nan]
    else:
        L2 = coord.stop_v2 - coord.start_v2

    # consider the first junction
    m1 = int(L1 / 3)
    if L1 % 3 == 0:
        spanning_id_range1 = set(range(max(m1 - k + 1, 0), m1))
    else:
        spanning_id_range1 = set(range(max(m1 - k + 1, 0), m1 + 1))

    # consider the second junction
    if coord.start_v3 is None:
        spanning_id_range2 = [np.nan]
        spanning_id_range1_2 = [np.nan]
    else:
        m2 = int((L1 + L2) / 3)
        if (L1 + L2) % 3 == 0:
            spanning_id_range2 = set(range(max(m2 - k + 1, 0), m2))
        else:
            spanning_id_range2 = set(range(max(m2 - k + 1, 0), m2 + 1))
        if (spanning_id_range1 and spanning_id_range2):
            spanning_id_range1_2 = spanning_id_range1.intersection(spanning_id_range2)
        else:
            spanning_id_range1_2 = []

    return spanning_id_range1, spanning_id_range2, spanning_id_range1_2

def get_and_write_background_peptide_and_kmer(peptide_set, kmer_set, gene, ref_mut_seq, gene_table, countinfo, kmer_length, all_read_frames):
    """Calculate the peptide translated from the complete transcript instead of single exon pairs

    Parameters
    ----------
    gene: Object. Created by SplAdder
    ref_mut_seq: List(str). Reference sequence of certain chromosome.
    gene_table: Namedtuple GeneTable, store the gene-transcript-cds mapping tables derived
       from .gtf file. has attribute ['gene_to_cds_begin', 'ts_to_cds', 'gene_to_cds']
    countinfo: NamedTuple containing SplAdder counts
    Idx: Namedtuple containing sample and gene index information
    kmer_length: int. length of kmer requested

    Returns
    -------
    background_peptide_list: List[str]. List of all the peptide translated from the given
       splicegraph and annotation.
    """

    gene_to_transcript_table, transcript_cds_table = gene_table.gene_to_ts, gene_table.ts_to_cds
    # Generate a background peptide for every variant transcript
    for ts in gene_to_transcript_table[gene.name]:
        # No CDS entries for transcript in annotation file...
        if ts not in transcript_cds_table:
            continue
        cds_expr_list, cds_string, cds_peptide_list = get_full_peptide(gene, ref_mut_seq['background'],
                                                                       transcript_cds_table[ts], countinfo,
                                                                       mode='back', all_read_frames=all_read_frames)

        for cds_peptide in cds_peptide_list: #all_read_frames modes outputs several peptides when encountering a stop codon
            peptide = OutputPeptide(output_id=ts,
                              peptide=cds_peptide,
                              exons_coor=None,
                              strand=None,
                              read_frame_annotated=None)
            peptide_set.add(namedtuple_to_str(peptide, sep = '\t'))

            if len(cds_peptide) >= kmer_length:
                for j in range(len(cds_peptide) - kmer_length + 1):
                    kmer_set.add(cds_peptide[j:j + kmer_length])


def retrieve_kmer_coordinates(start_pos_kmer, k, strand, spanning_index1, spanning_index2, spanning_index1_2,
                              translation_sorted_coord):
    '''
    :param start_pos_kmer: int. position inside the peptide where the kmer starts
    :param k: int. len of the kmer
    :param spanning_index1: set. set of starting positions for which the kmer will cross junction 1
    :param spanning_index2: set. set of starting positions for which the kmer will cross junction 2
    :param spanning_index1_2: set. set of starting positions for which the kmer will cross junctions 1 and 2
    :param translation_sorted_coord: np.array. exon coordinates modified by reading frame in translation order
    :return: Coord namedtuples with:
    start_v1 = int. start genomic coordinate in exon 1 of the kmer (may not be exon 1 of the peptide)
    stop_v1 = int. stop genomic coordinate in exon 1 of the kmer (may not be exon 1 of the peptide)
    start_v2 = int. start genomic coordinate in exon 2 of the kmer (may not be exon 2 of the peptide)
    stop_v2 = int. stop genomic coordinate in exon 2 of the kmer (may not be exon 2 of the peptide)
    start_v3 = int. start genomic coordinate in exon 3 of the kmer
    stop_v3 = int. stop genomic coordinate in exon 3 of the kmer
    '''

    def shift_to_start(start_pos_kmer, strand, overhang=0, end_kmer=False):
        if end_kmer:
            shift = ((start_pos_kmer + k) * 3) + overhang
        else:
            shift = (start_pos_kmer * 3) + overhang
        if strand == '-':
            shift = -shift
        return shift

    start_E1_modi, stop_E1_modi, start_E2_modi, stop_E2_modi, start_E3_modi, stop_E3_modi = translation_sorted_coord

    # get the position(s) right after the junction(s)
    if spanning_index1:
        pos_after_jx1 = max(spanning_index1) + 1
        first_cross_jx1 = min(spanning_index1)
    else:
        pos_after_jx1 = np.nan
        first_cross_jx1 = np.nan
    if spanning_index2:
        pos_after_jx2 = max(spanning_index2) + 1
    else:
        pos_after_jx2 = np.nan

    # get the overhang
    len_E1_modi = abs(stop_E1_modi - start_E1_modi)
    len_E2_modi = abs(stop_E2_modi - start_E2_modi)
    if (len_E1_modi % 3):
        overhang1 = 3 - (len_E1_modi % 3)
    else:
        overhang1 = 0
    if ((len_E1_modi + len_E2_modi) % 3):
        overhang2 = 3 - ((len_E1_modi + len_E2_modi) % 3)
    else:
        overhang2 = 0

    # get the number of nucleotides from the kmer relative to the junction(s)
    nt_kmer_left_jx1 = ((pos_after_jx1 - start_pos_kmer) * 3 ) - overhang1 #number of nucleotides before jx1 from the kmer
    nt_kmer_right_jx1 = ((k - (pos_after_jx1 - start_pos_kmer)) * 3 ) + overhang1
    nt_kmer_left_jx2 = ((pos_after_jx2 - start_pos_kmer) * 3) - overhang2
    nt_kmer_right_jx2 = ((k - (pos_after_jx2 - start_pos_kmer)) * 3) + overhang2
    if strand == '-':
        nt_kmer_left_jx1, nt_kmer_right_jx1 = -nt_kmer_left_jx1, -nt_kmer_right_jx1
        nt_kmer_left_jx2, nt_kmer_right_jx2 = -nt_kmer_left_jx2, -nt_kmer_right_jx2



    start_kmer_E2, stop_kmer_E2, start_kmer_E3, stop_kmer_E3 = np.nan, np.nan, None, None #nan and None for function get_segment_expr
    # Crosses 2 junctions
    if start_pos_kmer in spanning_index1_2:
        start_kmer_E1 = stop_E1_modi - nt_kmer_left_jx1
        stop_kmer_E1 = stop_E1_modi
        start_kmer_E2 = start_E2_modi
        stop_kmer_E2 = stop_E2_modi
        start_kmer_E3 = start_E3_modi
        stop_kmer_E3 = start_E3_modi + nt_kmer_right_jx2

    # Crosses first junction only
    elif start_pos_kmer in spanning_index1:
        start_kmer_E1 = stop_E1_modi - nt_kmer_left_jx1
        stop_kmer_E1 = stop_E1_modi
        start_kmer_E2 = start_E2_modi
        stop_kmer_E2 = start_E2_modi + nt_kmer_right_jx1

    # Crosses second junction only
    elif start_pos_kmer in spanning_index2:
        start_kmer_E1 = stop_E2_modi - nt_kmer_left_jx2
        stop_kmer_E1 = stop_E2_modi
        start_kmer_E2 = start_E3_modi
        stop_kmer_E2 = start_E3_modi + nt_kmer_right_jx2

    # Not cross junction: Before first junction
    elif start_pos_kmer < first_cross_jx1 or np.isnan(pos_after_jx1):
        start_kmer_E1 = start_E1_modi + shift_to_start(start_pos_kmer, strand)
        stop_kmer_E1 = start_E1_modi + shift_to_start(start_pos_kmer, strand, end_kmer=True)

    # Not cross junction: After second junction
    elif start_pos_kmer >= pos_after_jx2:
        start_kmer_E1 = start_E3_modi + shift_to_start(start_pos_kmer - pos_after_jx2, strand, overhang2)
        stop_kmer_E1 = start_E3_modi + shift_to_start(start_pos_kmer - pos_after_jx2, strand, overhang2, end_kmer=True)

    # Not cross junction: After first junction and potentially before second junction
    else:
        start_kmer_E1 = start_E2_modi + shift_to_start(start_pos_kmer - pos_after_jx1, strand, overhang1)
        stop_kmer_E1 = start_E2_modi + shift_to_start(start_pos_kmer - pos_after_jx1, strand, overhang1, end_kmer=True)

    if strand == '+':
        kmer_coord = Coord(start_kmer_E1, stop_kmer_E1,
                           start_kmer_E2, stop_kmer_E2,
                           start_kmer_E3, stop_kmer_E3)
    else:
        kmer_coord = Coord(stop_kmer_E1, start_kmer_E1,
                           stop_kmer_E2, start_kmer_E2,
                           stop_kmer_E3, start_kmer_E3)

    return kmer_coord


def create_kmer(output_peptide, k, gene_kmer_coord, kmer_database=None):

    """Calculate the output kmer and the corresponding expression based on output peptide #TODO update docstring!!!

    Parameters
    ----------
    output_peptide: OutputPeptide. Filtered output_peptide_list.
    k: int. Specify k-mer length
    graph_output_samples_ids: samples list found in graph
    gene_kmer_coord: dictionary of {coord: kmers} for the gene
    kmer_database: Set of kmers to be removed on the fly,
    usually from a public database as uniprot. I and L equivalence is applied

    Returns
    -------
    updates the kmer_matrix
    """

    # Get the peptide coordinates in translation order (once for each peptide)
    sort_coord = np.array([item for item in list(output_peptide.exons_coor) if item is not None])
    if output_peptide.strand == '-':
        sort_coord = - np.sort(-sort_coord)
    else:
        sort_coord = np.sort(sort_coord)
    sort_coord = np.pad(sort_coord, (0, 6 - len(sort_coord)), 'constant', constant_values=(-9999)) # add padding to get 6 entries

    # Get the peptide positions which span the respective junctions
    spanning_index1, spanning_index2, spanning_index1_2 = get_spanning_index(output_peptide.exons_coor, k)

    if len(output_peptide.peptide) >= k:
        for j in range(len(output_peptide.peptide) - k + 1):
            kmer_peptide = output_peptide.peptide[j:j+k]

            check_database = ((not kmer_database) or (replace_I_with_L(kmer_peptide) not in kmer_database)) # remove on the fly peptides from a database
            
            if check_database:
                kmer_coord = retrieve_kmer_coordinates(j, k, output_peptide.strand, spanning_index1, spanning_index2,
                                                       spanning_index1_2, sort_coord)
                gene_kmer_coord.add((kmer_coord, kmer_peptide, output_peptide.read_frame_annotated))


def prepare_output_kmer(gene, idx, countinfo, seg_counts, edge_idxs, edge_counts,
                                     gene_kmer_coord, gene_annot_jx,
                                     graph_output_samples_ids,
                                     graph_samples, filepointer, out_dir, verbose=False):

    '''
    :param gene:
    :param idx:
    :param k:
    :param countinfo:
    :param seg_counts:
    :param edge_idxs:
    :param edge_counts:
    :param gene_kmer_coord:
    :param gene_annot_jx:
    :param graph_output_samples_ids:
    :param graph_samples:
    :param filepointer:
    :param out_dir:
    :param verbose:
    :return:
    '''
    kmer_matrix_edge = []
    kmer_matrix_segm = []
    for kmer_coord, kmer_peptide, rf_annot in gene_kmer_coord:
        k = len(kmer_peptide)

        # segment expression
        _, pos_expr_segm = get_segment_expr(gene, kmer_coord, countinfo, idx, seg_counts)
        sublist_seg = np.round(np.atleast_2d(pos_expr_segm[:, 0]).dot(pos_expr_segm[:, 1:]) / (k * 3), 2)
        sublist_seg = sublist_seg[0].tolist()

        # junction expression
        if (countinfo is not None) and (kmer_coord.start_v2 is not np.nan): # kmer crosses only one exon
            _, edges_expr = search_edge_metadata_segmentgraph(gene, kmer_coord, edge_idxs, edge_counts)

            sublist_jun = np.nanmin(edges_expr, axis=0) # always apply. The min has no effect if one junction only
            if graph_output_samples_ids is not None:
                sublist_jun = sublist_jun[graph_output_samples_ids]
            sublist_jun = sublist_jun.tolist()
        else:
            sublist_jun = []

        # Flags
        if kmer_coord.start_v2 is np.nan: # kmer crosses only one exon
            is_in_junction = False
            junction_annotated = False
        else: # kmer crosses at least one junction
            is_in_junction = True
            jx1 = ':'.join([ str(i) for i in np.sort(kmer_coord[:4])[1:3]])
            junction_annotated = jx1 in gene_annot_jx
            if kmer_coord.start_v3 is not None:
                jx2 = ':'.join([str(i) for i in np.sort(kmer_coord[2:])[1:3]])
                junction_annotated = (jx1 in gene_annot_jx) or (jx2 in gene_annot_jx)

        # create output data
        row_metadata = [kmer_peptide, ':'.join([str(coord) for coord in kmer_coord]),
                        is_in_junction, junction_annotated, rf_annot]
        if is_in_junction:
            kmer_matrix_edge.append(row_metadata + sublist_jun)
        else:
            kmer_matrix_segm.append(row_metadata + sublist_seg)

        # save output data per batch
        if len(kmer_matrix_segm) > 1000:
            save_kmer_matrix(None, kmer_matrix_segm, graph_samples, filepointer, out_dir, verbose=False)
            kmer_matrix_segm.clear()
        if len(kmer_matrix_edge) > 1000:
            save_kmer_matrix(kmer_matrix_edge, None, graph_samples, filepointer, out_dir, verbose)
            kmer_matrix_edge.clear()

    save_kmer_matrix(kmer_matrix_edge, kmer_matrix_segm, graph_samples, filepointer, out_dir, verbose=False)
