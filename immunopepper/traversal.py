"""Contains all the output computation based on gene splicegraph"""

from collections import defaultdict
import numpy as np

from immunopepper.filter import add_peptide_properties
from immunopepper.filter import add_kmer_properties
from immunopepper.filter import add_kmers
from immunopepper.filter import filter_redundant_junctions
from immunopepper.filter import junction_is_annotated
from immunopepper.filter import is_intron_in_junction_list
from immunopepper.filter import junction_tuple_is_annotated
from immunopepper.io_ import save_kmer_matrix
from immunopepper.mutations import get_mutated_sequence
from immunopepper.mutations import exon_to_mutations
from immunopepper.mutations import get_mut_comb
from immunopepper.mutations import exon_to_expression
from immunopepper.namedtuples import Coord
from immunopepper.namedtuples import OutputBackground
from immunopepper.namedtuples import OutputJuncPeptide
from immunopepper.namedtuples import OutputKmer
from immunopepper.namedtuples import OutputMetadata
from immunopepper.namedtuples import VertexPair
from immunopepper.preprocess import search_edge_metadata_segmentgraph
from immunopepper.translate import get_exhaustive_reading_frames
from immunopepper.translate import get_full_peptide
from immunopepper.translate import isolated_peptide_result
from immunopepper.translate import get_peptide_result
from immunopepper.translate import cross_peptide_result
from immunopepper.utils import get_segment_expr


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



def collect_vertex_pairs(gene=None, gene_info=None, ref_seq_file=None, chrm=None, idx=None, mutation=None, all_read_frames=False, disable_concat=False, kmer=None, filter_redundant=False):
    """Calculate the output peptide for every exon-pair in the splicegraph

       Parameters
       ----------
       gene: Object, returned by SplAdder.
       ref_seq: Str, reference sequnce of specific chromosome
       idx: Namedtuple Idx, has attribute idx.gene and idx.sample
       mutation: Namedtuple Mutation, store the mutation information of specific chromosome and sample.
           has the attribute ['mode', 'maf_dict', 'vcf_dict']
       disable_concat: bool, flag indicating whether to disable the concanation of vertices into triples
        kmer: bool, flag indicating whether to extract kmers from the current parse
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
                        peptide, modi_coord, flag, next_reading_frame = cross_peptide_result(read_frame_tuple, gene.strand, variant_comb, mutation.somatic_dict, ref_mut_seq, sg.vertices[:, prop_vertex], min_pos, all_read_frames)
                        orig_coord = Coord(sg.vertices[0, v_id], sg.vertices[1, v_id], sg.vertices[0, prop_vertex], sg.vertices[1, prop_vertex])
                        if (not flag.has_stop) and (not all_read_frames): #no propagation needed in all reading frame mode
                            reading_frame_dict[prop_vertex].add(next_reading_frame)
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
        for kmer_length in kmer: #
            concat_vertex_pair_list['3-exons_{}-mer'.format(kmer_length)] = collect_vertex_triples(gene, vertex_pair_list, kmer_length)
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


def get_and_write_peptide_and_kmer(peptide_dict=None, kmer_dict=None,
                         gene=None, all_vertex_pairs=None, ref_mut_seq=None, idx=None,
                         exon_som_dict=None, countinfo=None,
                         edge_idxs=None, edge_counts=None, seg_counts=None,
                         mutation=None, mut_count_id=None, table=None,
                         size_factor=None, junction_list=None,
                         filepointer=None,
                         output_silence=False, kmer=None,
                         cross_graph_expr=None, all_read_frames=None, graph_output_samples_ids=None,
                         graph_samples=None,out_dir=None, verbose_save=None):
    """

    Parameters
    ----------
    peptide_dict: Dict. (peptide sequence)|-> metadata
    kmer_dict: Dict. (kmer sequence)|-> metadata
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
    filepointer: namedtuple, contains the columns and paths of each file of interest
    output_silence: bool, flag indicating whether not to silence annotated peptides
    kmer: bool, flag indicating whether to output kmers for this parse
    out_dir: str, base direactory used for temporary files
    cross_graph_expr: bool, whether to generate the expression kmer matrix with all samples from graph
    graph_samples: list, samples contained in the splicing graph object
    """
    # check whether the junction (specific combination of vertices) also is annotated
    # as a  junction of a protein coding transcript
    junction_flag = junction_is_annotated(gene, table.gene_to_ts, table.ts_to_cds)
    som_exp_dict = exon_to_expression(gene, list(mutation.somatic_dict.keys()), countinfo, seg_counts, mut_count_id)
    kmer_matrix = [{}, {},{}]  # in cross sample mode, will contain unique kmers per gene (1), is_junction (2), segments expr per sample (3), junction expr per sample (4)


    ### iterate over all vertex pairs and translate
    for kmer_type, vertex_pairs in all_vertex_pairs.items():
        for ii,vertex_pair in enumerate(vertex_pairs):
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
                    if not peptide.mut[pep_idx] or not (peptide.mut[pep_idx] not in peptide.ref or mutation.mode == 'ref' or output_silence):
                        continue

                    new_output_id = ':'.join([gene.name, '_'.join([str(v) for v in vertex_list]), str(variant_id), str(tran_start_pos), kmer_type])
                    vertex_tuple_anno_flag = junction_tuple_is_annotated(junction_flag, vertex_list)
                    is_intron_in_junction_list_flag = is_intron_in_junction_list(gene.splicegraph, vertex_list, gene.strand, junction_list)

                    if variant_comb is not np.nan and som_exp_dict is not None:  # which means there exist mutations
                        seg_exp_variant_comb = [int(som_exp_dict[ipos]) for ipos in variant_comb]
                    else:
                        seg_exp_variant_comb = np.nan  # if no mutation or no count file,  the segment expression is .

                    # collect expression data
                    if  countinfo:
                        segment_expr, expr_list = get_segment_expr(gene, modi_coord, countinfo, idx, seg_counts, cross_graph_expr=cross_graph_expr)
                    else:
                        segment_expr, expr_list = np.nan, None
                    if countinfo and not flag.is_isolated and edge_counts is not None: ## Will flag is isolated overlap with edge_counts is None?
                        edge_expr = search_edge_metadata_segmentgraph(gene, modi_coord, edge_idxs, edge_counts, cross_graph_expr=cross_graph_expr)
                    else:
                        edge_expr = np.nan

                    ### Peptides
                    add_peptide_properties(peptide_dict, [OutputMetadata(peptide=peptide.mut[pep_idx],
                                       output_id=new_output_id,
                                       read_frame=vertex_pair.read_frame.read_phase,
                                       gene_name=gene.name,
                                       gene_chr=gene.chr,
                                       gene_strand=gene.strand,
                                       mutation_mode=mutation.mode,
                                       junction_annotated=vertex_tuple_anno_flag,
                                       has_stop_codon=int(flag.has_stop),
                                       is_in_junction_list=is_intron_in_junction_list_flag,
                                       is_isolated=int(flag.is_isolated),
                                       variant_comb=variant_comb,
                                       variant_seg_expr=seg_exp_variant_comb,
                                       modified_exons_coord=modi_coord,
                                       original_exons_coord=vertex_pair.original_exons_coord,
                                       vertex_idx=vertex_list,
                                       junction_expr=edge_expr,
                                       segment_expr=segment_expr,
                                       kmer_type=kmer_type
                                       )], skip_expr=cross_graph_expr)
                    variant_id += 1
                    output_peptide = OutputJuncPeptide(output_id= new_output_id,
                                                    peptide=peptide.mut[pep_idx],
                                                    exons_coor=modi_coord,
                                                    junction_expr=edge_expr)

                    ### kmers
                    if cross_graph_expr: #generate kmer x sample expression matrix for all samples in graph
                        kmer_matrix = create_output_kmer_cross_samples(output_peptide, kmer[0], expr_list, graph_output_samples_ids, kmer_matrix) # Only one kmer lengthsupported for this mode

                    else:
                        if kmer:
                            if '2-exons' in kmer_type: #generate sample kmers for each vertex pair and each kmer_length
                                for kmer_length in kmer:
                                    add_kmer_properties(kmer_dict[kmer_length],
                                                             create_output_kmer(output_peptide, kmer_length, expr_list)) #TODO create other type of kmer dictionnary, + other saving function for batch // OR SAVE STRAIGHT away # adapt kmer to list of expressins -- keep cros junction info

                            else: #generate sample kmers for each vertex triplets, only for the kmer_lengths that require it
                                kmer_length = int(kmer_type.split('_')[-1].split('-')[0])
                                add_kmer_properties(kmer_dict[kmer_length],
                                                         create_output_kmer(output_peptide, kmer_length, expr_list))

        if not gene.splicegraph.edges is None:
            gene.to_sparse()
    if cross_graph_expr:
        save_kmer_matrix(kmer_matrix, graph_samples, filepointer, compression=None, out_dir=out_dir, verbose=verbose_save)



def get_spanning_index(coord, k):
    L1 = coord.stop_v1-coord.start_v1
    if coord.start_v2 is np.nan:
        return [np.nan], [np.nan], [np.nan]
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
        spanning_id_range2 = [np.nan]
        spanning_id_range1_2 = [np.nan]
    else:
        m2 = int((L1 + L2) / 3)
        if (L1 + L2) % 3 == 0:
            spanning_id_range2 = range(max(m2 - k + 1, 0), m2)
        else:
            spanning_id_range2 = range(max(m2 - k + 1, 0), m2 + 1)
        if (spanning_id_range1 and spanning_id_range2):
            spanning_id_range1_2 = [shared for shared in range(spanning_id_range2[0], spanning_id_range1[-1] + 1)
                                    if (shared + k - 1 >= m2 + 1)]
        else:
            spanning_id_range1_2 = []

    return spanning_id_range1, spanning_id_range2, spanning_id_range1_2

def get_and_write_background_peptide_and_kmer(peptide_dict, kmer_dict, gene, ref_mut_seq, gene_table, countinfo, seg_counts, Idx, kmer, all_read_frames):
    """Calculate the peptide translated from the complete transcript instead of single exon pairs

    Parameters
    ----------
    gene: Object. Created by SplAdder
    ref_mut_seq: List(str). Reference sequence of certain chromosome.
    gene_table: Namedtuple GeneTable, store the gene-transcript-cds mapping tables derived
       from .gtf file. has attribute ['gene_to_cds_begin', 'ts_to_cds', 'gene_to_cds']
    countinfo: NamedTuple containing SplAdder counts
    Idx: Namedtuple containing sample and gene index information
    kmer: bool, flag indicating whether kmers from the parse should be output

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
        cds_expr_list, cds_string, cds_peptide_list = get_full_peptide(gene, ref_mut_seq['background'], transcript_cds_table[ts], countinfo, seg_counts, Idx, mode='back', all_read_frames=all_read_frames)

        for cds_peptide in cds_peptide_list: #all_read_frames modes outputs several peptides when encountering a stop codon
            peptide = OutputBackground(output_id=ts, peptide=cds_peptide)
            add_peptide_properties(peptide_dict, [peptide])
            if kmer:
                for kmer_length in kmer:
                    add_kmers(kmer_dict[kmer_length],
                              create_output_kmer(peptide, kmer_length, cds_expr_list))



def create_output_kmer(output_peptide, k, expr_list):
    """Calculate the output kmer and the corresponding expression based on output peptide

    Parameters
    ----------
    output_peptide: OutputJuncPeptide. Filtered output_peptide_list.
    k: int. Specify k-mer length
    expr_lists: List(Tuple). Filtered expr_list.

    Returns
    -------
    output_kmer_list: List(str). Each line is the output peptide and corresponding expression level.

    """
    output_kmer_list = []
    peptide = output_peptide.peptide
    peptide_head = output_peptide.output_id
    if hasattr(output_peptide,'exons_coor'):
        coord = output_peptide.exons_coor
        spanning_index1, spanning_index2, spanning_index1_2 = get_spanning_index(coord, k)
    else:
        spanning_index1, spanning_index2, spanning_index1_2 = [np.nan], [np.nan], [np.nan]
    if hasattr(output_peptide, 'junction_expr'):
        junction_count = output_peptide.junction_expr
    else:
        junction_count = np.nan
    # decide the kmer that spans over the cross junction
    if expr_list is None:
        expr_array = None
    else:
        expr_array = np.array([x[1] for x in expr_list for _ in range(int(x[0]))])
    if len(peptide) >= k:
        for j in range(len(peptide) - k + 1):
            kmer_peptide = peptide[j:j+k]
            if expr_array is None:
                kmer_peptide_expr = np.nan
            else:
                kmer_peptide_expr = np.round(np.mean(expr_array[j*3:(j+k)*3]), 2)
            if j in spanning_index1_2:
                is_in_junction = True
                kmer_junction_count = np.nanmin(junction_count)
            elif j in spanning_index1:
                is_in_junction = True
                kmer_junction_count = junction_count[0] if junction_count is not np.nan else np.nan
            elif j in spanning_index2 :
                is_in_junction = True
                kmer_junction_count = junction_count[1] if junction_count is not np.nan else np.nan
            else:
                is_in_junction = False
                kmer_junction_count = np.nan
            kmer = OutputKmer(kmer_peptide, peptide_head, kmer_peptide_expr, is_in_junction, kmer_junction_count)
            output_kmer_list.append(kmer)
    return output_kmer_list


def create_output_kmer_cross_samples(output_peptide, k, segm_expr_list, graph_output_samples_ids, kmer_matrix):
    """Calculate the output kmer and the corresponding expression based on output peptide

    Parameters
    ----------
    output_peptide: OutputJuncPeptide. Filtered output_peptide_list.
    k: int. Specify k-mer length
    segm_expr_list: List(Tuple). Filtered expr_list.
    graph_samples: samples list found in graph
    kmer_matrix: [will contain unique kmers per gene (1), is_junction (2), segments expr per sample (3), junction expr per sample (4)]


    Returns
    -------
    updates the kmer_matrix
    """

    peptide = output_peptide.peptide
    positions = np.cumsum(segm_expr_list[:, 0])

    if hasattr(output_peptide,'exons_coor'):
        coord = output_peptide.exons_coor
        spanning_index1, spanning_index2, spanning_index1_2 = get_spanning_index(coord, k)
    else:
        spanning_index1, spanning_index2, spanning_index1_2 = [np.nan], [np.nan], [np.nan]
    if hasattr(output_peptide, 'junction_expr'):
        junction_count = output_peptide.junction_expr
    else:
        junction_count = np.nan


    if len(peptide) >= k:
        for j in range(len(peptide) - k + 1):
            kmer_peptide = peptide[j:j+k]

            # junction expression
            is_in_junction = True
            if j in spanning_index1_2:
                sublist_jun = np.nanmin(np.array(junction_count), axis=0)[graph_output_samples_ids] if junction_count is not np.nan \
                    else np.empty((len(graph_output_samples_ids),)) * np.nan
            elif j in spanning_index1:
                sublist_jun = junction_count[0][graph_output_samples_ids] if junction_count is not np.nan \
                    else np.empty((len(graph_output_samples_ids),)) * np.nan
            elif j in spanning_index2:
                sublist_jun = junction_count[1][graph_output_samples_ids] if junction_count is not np.nan \
                    else np.empty((len(graph_output_samples_ids),)) * np.nan
            else:
                is_in_junction = False
                sublist_jun = np.empty((len(graph_output_samples_ids),)) * np.nan

            # segment expression
            if segm_expr_list is None:
                sublist_seg = [np.nan] * len(graph_output_samples_ids)
            else:
                W = np.zeros(segm_expr_list.shape[0])
                left = min(np.where(positions >= (j * 3))[0])
                right = min(np.where(positions >= (j + k) * 3)[0])
                if left != right:
                    W[left] = positions[left] - j * 3
                    for i in np.arange(left + 1, right, 1):
                        W[i] = segm_expr_list[i, 0]
                    W[right] = (j + k) * 3 - positions[right - 1]
                else:
                    W[left] = (j + k) * 3 - j * 3
                if j == 0:
                    sublist_seg =  np.atleast_2d(W[left:right+1]).dot(segm_expr_list[left:right + 1, 1:]) / (k*3)
                else:
                    delta = W - W_past
                    sublist_seg = sublist_seg
                    for pos_ in np.arange(left_past, right + 1, 1) :
                        sublist_seg += (delta[pos_] * segm_expr_list[pos_, 1:]) / (k*3)
                sublist_seg = sublist_seg.flatten().tolist()
                left_past = left
                W_past = W

            # update the cross samples matrix
            if sum(np.isnan(sublist_seg)) != len(sublist_seg) or sum(np.isnan(sublist_jun)) != len(sublist_jun):
                if kmer_peptide not in kmer_matrix[0]:
                    kmer_matrix[0][kmer_peptide] = is_in_junction
                    kmer_matrix[1][kmer_peptide] = np.round(sublist_seg, 2)
                    kmer_matrix[2][kmer_peptide] = sublist_jun
                else:
                    kmer_matrix[0][kmer_peptide] = max(kmer_matrix[0][kmer_peptide] ,is_in_junction )
                    kmer_matrix[1][kmer_peptide] = np.nanmax(np.array( [kmer_matrix[1][kmer_peptide], np.round(sublist_seg, 2)]), axis = 0)# make unique per gene with maximum
                    kmer_matrix[2][kmer_peptide] = np.nanmax(np.array( [kmer_matrix[2][kmer_peptide], sublist_jun]), axis = 0)


    return kmer_matrix





