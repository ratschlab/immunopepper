"""Contains all the output computation based on gene splicegraph"""

import logging
import h5py
import numpy as np

from .constant import NOT_EXIST
from .filter import get_filtered_metadata_list
from .filter import junction_is_annotated
from .filter import junction_is_in_given_list
from .filter import junction_tuple_is_annotated
from .filter import peptide_is_annotated
from .io import convert_namedtuple_to_str
from .io import write_namedtuple_list
from .mutations import apply_germline_mutation
from .mutations import apply_somatic_mutation
from .mutations import get_exon_som_dict
from .mutations import get_mut_comb
from .mutations import get_som_expr_dict
from .namedtuples import Coord
from .namedtuples import OutputBackground
from .namedtuples import OutputJuncPeptide
from .namedtuples import OutputKmer
from .namedtuples import OutputMetadata
from .namedtuples import VertexPair
from .preprocess import search_edge_metadata_segmentgraph
from .translate import get_full_peptide
from .translate import isolated_peptide_result
from .translate import get_peptide_result
from .translate import cross_peptide_result
from .utils import get_segment_expr


def collect_vertex_pairs(gene=None, gene_info=None, ref_seq_file=None, chrm=None, idx=None, mutation=None, disable_concat=False, kmer=None, filter_redundant=False):
    """Calculte the output peptide for every exon-pairs in the splicegraph

       Parameters
       ----------
       gene: Object, returned by SplAdder.
       ref_seq: Str, reference sequnce of specific chromosome
       idx: Namedtuple Idx, has attribute idx.gene and idx.sample
       mutation: Namedtuple Mutation, store the mutation information of specific chromosome and sample.
           has the attribute ['mode', 'maf_dict', 'vcf_dict']

       Returns
       -------
       vertex_pairs: List of VertexPair.
       ref_mut_seq: Dict. (sequence_type) -> list[char].
       exon_som_dict: Dict. (exon_id) |-> (mutation_postion)

       """

    gene.from_sparse()
    sg = gene.splicegraph
    min_pos = np.min(sg.vertices[0])
    max_pos = np.max(sg.vertices[1])

    output_id = 0

    # apply germline mutation
    # when germline mutation is applied, background_seq != ref_seq
    # otherwise, background_seq = ref_seq
    ref_mut_seq = apply_germline_mutation(ref_sequence_file=ref_seq_file,
                                          chrm=chrm,
                                          pos_start=min_pos,
                                          pos_end=max_pos,
                                          mutation_sub_dict=mutation.germline_mutation_dict)

    # apply somatic mutation
    # exon_som_dict: (exon_id) |-> (mutation_postion)
    exon_som_dict = None
    if mutation.somatic_mutation_dict is not None:
        exon_som_dict = get_exon_som_dict(gene, mutation.somatic_mutation_dict)

    vertex_pair_list = []
    reading_frame_dict = dict(gene_info.reading_frames)

    for v_id in gene_info.vertex_order:
        n_read_frames = len(reading_frame_dict[v_id])
        if n_read_frames == 0:  # no cds start, skip the vertex
            continue
        if len(gene_info.vertex_succ_list[v_id]) == 0: # if no successive vertex, we add a flag NOT_EXIST, translate and output it
            gene_info.vertex_succ_list[v_id].append(NOT_EXIST)

        for prop_vertex in gene_info.vertex_succ_list[v_id]:
            vertex_list = [v_id, prop_vertex]
            mut_seq_comb = get_mut_comb(exon_som_dict,vertex_list)
            for read_frame_tuple in sorted(reading_frame_dict[v_id]):
                has_stop_flag = True
                for variant_comb in mut_seq_comb:  # go through each variant combination
                    logging.debug(' '.join([str(v_id), str(prop_vertex), str(variant_comb), str(read_frame_tuple.read_phase)]))
                    if prop_vertex != NOT_EXIST:
                        peptide, modi_coord, flag, next_reading_frame = cross_peptide_result(read_frame_tuple, gene.strand, variant_comb, mutation.somatic_mutation_dict, ref_mut_seq, sg.vertices[:, prop_vertex], min_pos)
                        orig_coord = Coord(sg.vertices[0, v_id], sg.vertices[1, v_id], sg.vertices[0, prop_vertex], sg.vertices[1, prop_vertex])
                        if not flag.has_stop:
                            reading_frame_dict[prop_vertex].add(next_reading_frame)
                    else:
                        peptide, modi_coord, flag = isolated_peptide_result(read_frame_tuple, gene.strand, variant_comb, mutation.somatic_mutation_dict, ref_mut_seq, min_pos)
                        orig_coord = Coord(sg.vertices[0, v_id],sg.vertices[1, v_id], NOT_EXIST, NOT_EXIST)
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

    if  disable_concat:
        concat_vertex_pair_list = []
    else:
        concat_vertex_pair_list = collect_vertex_triples(gene, vertex_pair_list, kmer)
    if filter_redundant:
        vertex_pair_list = get_filtered_metadata_list(vertex_pair_list, gene.strand)
    vertex_pair_list += concat_vertex_pair_list

    return vertex_pair_list, ref_mut_seq, exon_som_dict


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
                        and vp.vertex_idxs[1] != NOT_EXIST]
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


def get_and_write_peptide_and_kmer(gene=None, vertex_pairs=None, background_pep_list=None, ref_mut_seq=None, idx=None,
                         exon_som_dict=None, countinfo=None, mutation=None,table=None,
                         size_factor=None, junction_list=None, output_silence=False, kmer=None):
    """

    Parameters
    ----------
    gene: Object, returned by SplAdder.
    vertex_pairs: List of VertexPair
    background_pep_list: List[str]. List of all the peptide translated from the given splicegraph and annotation.
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

    """
    # check whether the junction (specific combination of vertices) also is annotated
    # as a  junction of a protein coding transcript
    junction_flag = junction_is_annotated(gene, table.gene_to_ts, table.ts_to_cds)
    som_exp_dict = get_som_expr_dict(gene, list(mutation.somatic_mutation_dict.keys()), countinfo, idx)

    ### collect the relevant count infor for the current gene
    if countinfo:
        gidx = countinfo.gene_idx_dict[gene.name]
        edge_gene_idxs = np.arange(countinfo.gene_id_to_edgerange[gidx][0], countinfo.gene_id_to_edgerange[gidx][1])
        with h5py.File(countinfo.h5fname, 'r') as h5f:
            edge_idxs = h5f['edge_idx'][edge_gene_idxs].astype('int')
            edge_counts = h5f['edges'][edge_gene_idxs, idx.sample] 
            seg_gene_idxs = np.arange(countinfo.gene_id_to_segrange[gidx][0], countinfo.gene_id_to_segrange[gidx][1])
            seg_counts = h5f['segments'][seg_gene_idxs, idx.sample] 

    output_metadata_list = []
    output_peptide_list = []
    output_kmer_lists = []

    ### iterate over all vertex pairs and translate
    for ii,vertex_pair in enumerate(vertex_pairs):
        modi_coord = vertex_pair.modified_exons_coord
        vertex_list = vertex_pair.vertex_idxs
        tran_start_pos = modi_coord.start_v1 if gene.strand == '+' else modi_coord.stop_v1
        mut_seq_comb = get_mut_comb(exon_som_dict, vertex_pair.vertex_idxs)

        variant_id = 0
        for variant_comb in mut_seq_comb:  # go through each variant combination
            peptide,flag = get_peptide_result(vertex_pair, gene.strand, variant_comb, mutation.somatic_mutation_dict, ref_mut_seq, np.min(gene.splicegraph.vertices))

            # If cross junction peptide has a stop-codon in it, the frame
            # will not be propagated because the read is truncated before it reaches the end of the exon.
            # also in mutation mode, only output the case where ref is different from mutated
            if not peptide.mut or not (peptide.mut != peptide.ref or mutation.mode == 'ref' or output_silence):
                continue

            new_output_id = ':'.join([gene.name, '_'.join([str(v) for v in vertex_list]), str(variant_id), str(tran_start_pos)])
            peptide_is_annotated_flag = len(peptide_is_annotated(background_pep_list, peptide.mut))
            vertex_tuple_anno_flag = junction_tuple_is_annotated(junction_flag, vertex_list)
            junction_is_in_given_list_flag = junction_is_in_given_list(gene.splicegraph, vertex_list, gene.strand, junction_list)

            if variant_comb != NOT_EXIST and som_exp_dict is not None:  # which means there exist mutations
                seg_exp_variant_comb = [int(som_exp_dict[ipos]) for ipos in variant_comb]
            else:
                seg_exp_variant_comb = NOT_EXIST  # if no mutation or no count file,  the segment expression is .

            # collect expression data
            if  countinfo:
                segment_expr, expr_list = get_segment_expr(gene, modi_coord, countinfo, idx, seg_counts)
            else:
                segment_expr, expr_list = NOT_EXIST, None
            if countinfo and not flag.is_isolated:
                edge_expr = search_edge_metadata_segmentgraph(gene, modi_coord, countinfo, idx, edge_idxs, edge_counts)
            else:
                edge_expr = NOT_EXIST

            output_metadata_list.append(OutputMetadata(output_id=new_output_id,
                                                       read_frame=vertex_pair.read_frame.read_phase,
                                                       gene_name=gene.name,
                                                       gene_chr=gene.chr,
                                                       gene_strand=gene.strand,
                                                       mutation_mode=mutation.mode,
                                                       peptide_annotated=peptide_is_annotated_flag,
                                                       junction_annotated=vertex_tuple_anno_flag,
                                                       has_stop_codon=int(flag.has_stop),
                                                       is_in_junction_list=junction_is_in_given_list_flag,
                                                       is_isolated=int(flag.is_isolated),
                                                       variant_comb=variant_comb,
                                                       variant_seg_expr=seg_exp_variant_comb,
                                                       modified_exons_coord=modi_coord,
                                                       original_exons_coord=vertex_pair.original_exons_coord,
                                                       vertex_idx=vertex_list,
                                                       junction_expr=edge_expr,
                                                       segment_expr=segment_expr
            ))
            variant_id += 1
            output_peptide = OutputJuncPeptide(output_id='>' + new_output_id,
                                            id=new_output_id,
                                            peptide=peptide.mut,
                                            exons_coor=modi_coord,
                                            junction_count=edge_expr)
            output_peptide_list.append(output_peptide)
            if kmer > 0:
                output_kmer_lists.append(create_output_kmer(output_peptide, kmer, expr_list))

    if not gene.splicegraph.edges is None:
        gene.to_sparse()
    gene.processed = True

    return output_metadata_list, output_peptide_list, output_kmer_lists


def get_and_write_background_peptide_and_kmer(gene, ref_mut_seq, gene_table, countinfo, Idx, kmer):
    """Calculate the peptide translated from the complete transcript instead of single exon pairs

    Parameters
    ----------
    gene: Object. Created by SplAdder
    ref_mut_seq: List(str). Reference sequence of certain chromosome.
    gene_table: Namedtuple GeneTable, store the gene-transcript-cds mapping tables derived
       from .gtf file. has attribute ['gene_to_cds_begin', 'ts_to_cds', 'gene_to_cds']
    countinfo: NamedTuple containing SplAdder counts

    Returns
    -------
    peptide_list: List[str]. List of all the peptide translated from the given
       splicegraph and annotation.
    (ts_list): List[str]. List of all the transcript indicated by the  annotation file
        can be used to generate artifical reads.
    """
    gene_to_transcript_table, transcript_cds_table = gene_table.gene_to_ts, gene_table.ts_to_cds
    background_peptide_list = []
    output_kmer_lists = []
    # Generate a background peptide for every variant transcript
    for ts in gene_to_transcript_table[gene.name]:
        # No CDS entries for transcript in annotation file...
        if ts not in transcript_cds_table:
            continue
        cds_expr_list, cds_string, cds_peptide = get_full_peptide(gene, ref_mut_seq['background'], transcript_cds_table[ts], countinfo, Idx, mode='back')
        peptide = OutputBackground(ts, cds_peptide)
        background_peptide_list.append(peptide)
        if kmer > 0:
            output_kmer_lists.append(create_output_kmer(peptide, kmer, cds_expr_list))
    gene.processed = True
    return background_peptide_list, output_kmer_lists


def create_output_kmer(output_peptide, k, expr_list):
    """Calculate the output kmer and the corresponding expression based on output peptide

    Parameters
    ----------
    peptide_list: OutputJuncPeptide. Filtered output_peptide_list.
    expr_lists: List(Tuple). Filtered expr_list.
    k: int. Specify k-mer length

    Returns
    -------
    output_list: List(str). Each line is the output peptide and corresponding expression level.

    """
    def get_spanning_index(coord, k):
        L1 = coord.stop_v1-coord.start_v1
        if coord.start_v2 == NOT_EXIST:
            return [NOT_EXIST],[NOT_EXIST]
        else:
            L2 = coord.stop_v2 - coord.start_v2

        # consider the first junction
        m1 = int(L1 / 3)
        if L1 % 3 == 0:
            spanning_id_range1 = range(max(m1 - k + 1, 0), m1)
        else:
            spanning_id_range1 = range(max(m1 - k + 1, 0), m1 + 1)

        # consider the second junction
        if coord.start_v3 is None:
            spanning_id_range2 = [NOT_EXIST]
        else:
            m2 = int((L1 + L2) / 3)
            if (L1 + L2) % 3 == 0:
                spanning_id_range2 = range(max(m2 - k + 1, 0), m2)
            else:
                spanning_id_range2 = range(max(m2 - k + 1, 0), m2 + 1)

        return spanning_id_range1, spanning_id_range2

    output_kmer_list = []
    peptide = output_peptide.peptide
    peptide_head = output_peptide.id
    if hasattr(output_peptide,'exons_coor'):
        coord = output_peptide.exons_coor
        spanning_index1, spanning_index2 = get_spanning_index(coord, k)
    else:
        spanning_index1, spanning_index2 = [NOT_EXIST], [NOT_EXIST]
    if hasattr(output_peptide, 'junction_count'):
        junction_count = output_peptide.junction_count
    else:
        junction_count = NOT_EXIST
    # decide the kmer that spans over the cross junction
    if expr_list is None:
        expr_array = None
    else:
        expr_array = np.array([x[1] for x in expr_list for _ in range(int(x[0]))])
    if len(peptide) >= k:
        for j in range(len(peptide) - k + 1):
            kmer_peptide = peptide[j:j+k]
            if expr_array is None:
                kmer_peptide_expr = NOT_EXIST
            else:
                kmer_peptide_expr = np.round(np.mean(expr_array[j*3:(j+k)*3]), 2)
            if j in spanning_index1:
                is_in_junction = True
                kmer_junction_count = junction_count[0] if junction_count != NOT_EXIST else NOT_EXIST
            elif j in spanning_index2 :
                is_in_junction = True
                kmer_junction_count = junction_count[1] if junction_count != NOT_EXIST else NOT_EXIST
            else:
                is_in_junction = False
                kmer_junction_count = NOT_EXIST
            kmer = OutputKmer(kmer_peptide, peptide_head, kmer_peptide_expr, is_in_junction, kmer_junction_count)
            output_kmer_list.append(kmer)
    return output_kmer_list

