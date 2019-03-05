"""Contains all the output computation based on gene splicegraph"""
from __future__ import print_function

import numpy as np
import scipy as sp

from immuno_filter import junction_is_annotated, peptide_match, find_background_peptides
from immuno_mutation import apply_germline_mutation,get_exon_som_dict,get_som_expr_dict,get_mut_comb,apply_somatic_mutation
from utils import cross_peptide_result,is_isolated_cds,isolated_peptide_result,is_in_junction_list,get_segment_expr
from immuno_preprocess import search_edge_metadata_segmentgraph
from constant import NOT_EXIST
from immuno_nametuple import Output_metadata, Output_junc_peptide, Output_kmer

def calculate_output_peptide(gene=None, ref_seq=None, idx=None,
                      segments=None, edges=None, mutation=None,
                             table=None, debug=False, output_silence=False,
                             size_factor=None, junction_list=None):
    """Calculte the output peptide for every exon-pairs in the splicegraph
       Parameters
       ----------
       gene: Object, returned by SplAdder.
       ref_seq: Str, reference sequnce of specific chromosome
       idx: Namedtuple Idx, has attribute idx.gene and idx.sample
       segments: Namedtuple Segments, store segment expression information from count.hdf5.
           has attribute ['expr', 'lookup_table'].
       edges: Namedtuple Edges, store edges expression information from count.hdf5.
           has attribute ['expr','lookup_table']
       table: Namedtuple GeneTable, store the gene-transcript-cds mapping tables derived
           from .gtf file. has attribute ['gene_to_cds_begin', 'ts_to_cds', 'gene_to_cds']
       debug: bool. More detailed information will be printed when debugging.
       output_silence: bool. If set true, mutated peptide will be output even it is the same as referencd peptide
       size_factor: Scalar. To adjust the expression counts based on the external file `libsize.tsv`
       junction_list: List. Work as a filter to indicate some exon pair has certain
           ordinary intron which can be ignored further.
       mutation: Namedtuple Mutation, store the mutation information of specific chromosome and sample.
           has the attribute ['mode', 'maf_dict', 'vcf_dict']
       Returns
       -------
       output_peptide_list: List[Output_junc_peptide]. Contain all the possible output peptide in the given splicegraph.
       output_metadata_list: List[Output_metadata]. Contain the correpsonding medata data for each output peptide.
       output_background_pep_list: List[Output_background]. Contain the background peptide for each transcript.
       expr_lists: List[List(Tuple(int,float))]. Contain the segment expression data for each output peptide.
       back_expr_lists: List[List(Tuple(int,float))]. Contain the segment expression data for background output.
       total_expr: Float. The sum of all the expression counts which will be used for generating libsize.tsv
       """

    sg = gene.splicegraph
    gene.from_sparse()
    total_expr = 0

    output_id = 0
    output_peptide_list = []
    output_metadata_list = []
    expr_lists = []

    # apply germline mutation
    # when germline mutation is applied, background_seq != ref_seq
    # otherwise, background_seq = ref_seq
    pos_start = np.min(sg.vertices[0])
    pos_end = np.max(sg.vertices[1])
    ref_mut_seq = apply_germline_mutation(ref_sequence=ref_seq, pos_start=pos_start, pos_end=pos_end,
                                           mutation_sub_dic_vcf=mutation.vcf_dict)

    # apply somatic mutation
    # som_exp_dict: (mutation_position) |-> (expression)
    # exon_som_dict: (exon_id) |-> (mutation_postion)
    som_exp_dict, exon_som_dict = None,None
    if mutation.maf_dict is not None:
        exon_som_dict = get_exon_som_dict(gene, mutation.maf_dict.keys())
        if segments is not None:
            som_exp_dict = get_som_expr_dict(gene, mutation.maf_dict.keys(), segments, idx)

    # find background peptide
    # if no germline mutation is applies, germline key still exists, equals to reference.
    # return the list of the background peptide for each transcript
    background_pep_list, back_expr_lists = find_background_peptides(gene, ref_mut_seq['background'], table.gene_to_ts, table.ts_to_cds, segments, idx)

    # check whether the junction (specific combination of vertices) also is annotated
    # as a  junction of a protein coding transcript
    junction_flag = junction_is_annotated(gene, table.gene_to_ts, table.ts_to_cds)
    reading_frame_dict = dict(sg.reading_frames)
    for v_id in gene.vertex_order:
        n_read_frames = len(reading_frame_dict[v_id])
        if n_read_frames == 0:  # no cds start, skip the vertex
            continue
        if len(gene.vertex_succ_list[v_id]) == 0: # if no successive vertex, we add a flag NOT_EXIST, translate and output it
            gene.vertex_succ_list[v_id].append(NOT_EXIST)
        for prop_vertex in gene.vertex_succ_list[v_id]:
            mut_seq_comb = get_mut_comb(exon_som_dict, v_id, prop_vertex)
            for variant_comb in mut_seq_comb:  # go through each variant combination
                for read_frame_tuple in sorted(reading_frame_dict[v_id]):
                    if debug:
                        print(v_id, prop_vertex, variant_comb, read_frame_tuple)
                    peptide_weight = 1.0 / n_read_frames
                    if prop_vertex != NOT_EXIST:
                        peptide, coord, flag, next_reading_frame = cross_peptide_result(read_frame_tuple, gene.strand, variant_comb, mutation.maf_dict, ref_mut_seq, sg.vertices[:, prop_vertex])
                        if not flag.has_stop:
                            reading_frame_dict[prop_vertex].add(next_reading_frame)
                    else:
                        peptide, coord, flag = isolated_peptide_result(read_frame_tuple, gene.strand, variant_comb, mutation.maf_dict,ref_mut_seq)

                    # If cross junction peptide has a stop-codon in it, the frame
                    # will not be propagated because the read is truncated before it reaches the end of the exon.
                    # also in mutation mode, only output the case where ref is different from mutated
                    if peptide.mut != peptide.ref or mutation.mode == 'ref' or output_silence:
                        match_ts_list = peptide_match(background_pep_list, peptide.mut)
                        peptide_is_annotated = len(match_ts_list)
                        if not flag.is_isolated:
                            junction_anno_flag = int(junction_flag[v_id, prop_vertex])
                            if junction_list is not None:
                                if gene.strand == '+':
                                    junctionOI_flag = is_in_junction_list(sg.vertices[:, v_id], sg.vertices[:, prop_vertex], gene.strand, junction_list)
                                else:
                                    junctionOI_flag = is_in_junction_list(sg.vertices[:, prop_vertex], sg.vertices[:, v_id], gene.strand, junction_list)
                            else:
                                junctionOI_flag = NOT_EXIST
                        else:
                            junction_anno_flag = NOT_EXIST
                            junctionOI_flag = NOT_EXIST
                        # Write the variant gene into the FASTA FP together with the donor ID
                        str_variant_comb = [str(ipos) for ipos in variant_comb]

                        if variant_comb != NOT_EXIST and som_exp_dict is not None:  # which means there do exist some mutation
                            seg_exp_variant_comb = [str(int(som_exp_dict[ipos])) for ipos in variant_comb]
                        else:
                            seg_exp_variant_comb = NOT_EXIST  # if no mutation or no count file,  the segment expression is .

                        if segments is not None:
                            segment_expr, expr_list = get_segment_expr(gene, coord, segments, idx)
                        else:
                            segment_expr, expr_list = NOT_EXIST, NOT_EXIST

                        # create namedtuple
                        gene_outputid = str(idx.gene)+'.'+str(output_id)
                        # deal with expression data
                        if edges is not None and not flag.is_isolated:
                            sorted_pos = sp.sort(np.array([coord.start_v1, coord.stop_v1, coord.start_v2, coord.stop_v2]))
                            edge_expr = search_edge_metadata_segmentgraph(gene, sorted_pos, edges, idx)
                            total_expr += edge_expr
                            #edge_expr = edge_expr*size_factor
                        else:
                            edge_expr = NOT_EXIST
                        output_metadata = Output_metadata(output_id=gene_outputid,
                                                          read_frame=read_frame_tuple.read_phase,
                                                          gene_name=gene.name,
                                                          gene_chr=gene.chr,
                                                          gene_strand=gene.strand,
                                                          mutation_mode=mutation.mode,
                                                          peptide_weight="{:.3f}".format(peptide_weight),
                                                          peptide_annotated=peptide_is_annotated,
                                                          junction_annotated=junction_anno_flag,
                                                          has_stop_codon=int(flag.has_stop),
                                                          is_in_junction_list=junctionOI_flag,
                                                          is_isolated=int(flag.is_isolated),
                                                          variant_comb=str_variant_comb,
                                                          variant_seg_expr=seg_exp_variant_comb,
                                                          exons_coor=coord,
                                                          vertex_idx=str(v_id) + ',' + str(prop_vertex),
                                                          junction_expr=edge_expr,
                                                          segment_expr=segment_expr
                        )
                        output_metadata_list.append(output_metadata)
                        output_peptide = Output_junc_peptide(output_id='>'+str(idx.gene)+'.'+str(output_id),
                                                        id=gene.name+'_'+str(v_id)+'_'+str(prop_vertex),
                                                        peptide=peptide.mut,exons_coor=coord)
                        output_peptide_list.append(output_peptide)
                        expr_lists.append(expr_list)
                        output_id += 1
    if not sg.edges is None:
        gene.to_sparse()

    gene.processed = True
    return output_peptide_list,output_metadata_list,background_pep_list,expr_lists,back_expr_lists,total_expr


def create_output_kmer(peptide_list, expr_lists, k):
    """Calculate the output kmer and the corresponding expression based on output peptide

    Parameters
    ----------
    peptide_list: List(Output_junc_peptide). Filtered output_peptide_list.
    expr_lists: List(List(Tuple)). Filtered expr_list.
    k: int. Specify k-mer length

    Returns
    -------
    output_list: List(str). Each line is the output peptide and corresponding expression level.

    """
    def change_expr_lists_to_array(expr_list):
        array = []
        if NOT_EXIST in expr_list:
            array = [NOT_EXIST]
        else:
            for item in expr_list:
                length = item[0]
                expr = item[1]
                array.extend([expr]*length)
        return np.array(array)

    def get_spanning_index(coord, k):
        L1 = coord.stop_v1-coord.start_v1
        if coord.start_v2 == NOT_EXIST:
            spanning_id_range  = NOT_EXIST
            return spanning_id_range
        else:
            L2 = coord.stop_v2-coord.start_v2
        m = L1 / 3
        if L1%3 == 0:
            spanning_id_range = range(max(m-k+1,0),m)
        else:
            spanning_id_range = range(max(m-k+1,0),m+1)
        return spanning_id_range

    assert len(peptide_list) == len(expr_lists)
    output_list = []
    for i in range(len(peptide_list)):
        peptide = peptide_list[i].peptide
        peptide_head = peptide_list[i].id
        if hasattr(peptide_list[i],'exons_coor'):
            coord = peptide_list[i].exons_coor
            spanning_index = get_spanning_index(coord,k)
        else:
            spanning_index = NOT_EXIST
        # decide the kmer that spans over the cross junction
        expr_array = change_expr_lists_to_array(expr_lists[i])
        if len(peptide) >= k:
            for j in range(len(peptide)-k+1):
                kmer_peptide = peptide[j:j+k]
                if NOT_EXIST in expr_array:
                    kmer_peptide_expr = NOT_EXIST
                else:
                    kmer_peptide_expr = np.round(np.mean(expr_array[j*3:(j+k)*3]),2)
                if spanning_index is NOT_EXIST:
                    is_in_junction = NOT_EXIST
                else:
                    is_in_junction = j in spanning_index
                kmer = Output_kmer(kmer_peptide,peptide_head,kmer_peptide_expr,is_in_junction)
                output_list.append(kmer)
        else:
            kmer_peptide = peptide
            if spanning_index is NOT_EXIST:
                is_in_junction = False
            else:
                is_in_junction = True
            if NOT_EXIST in expr_array:
                kmer_peptide_expr = NOT_EXIST
            else:
                kmer_peptide_expr = np.round(np.mean(expr_array),2)
            kmer = Output_kmer(kmer_peptide, peptide_head, kmer_peptide_expr,is_in_junction)
            output_list.append(kmer)
    return output_list
