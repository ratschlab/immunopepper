from __future__ import print_function
# External libraries
import numpy as np
import scipy as sp
# immuno module
from immuno_filter import junction_is_annotated, peptide_match, find_background_peptides
from immuno_mutation import apply_germline_mutation,get_exon_som_dict,get_som_expr_dict,get_mut_comb
from utils import cross_peptide_result,is_isolated_cds,isolated_peptide_result,is_output_redundant,is_in_junction_list
from immuno_preprocess import search_edge_metadata_segmentgraph


# Optimized annotation code that does not loop over the annotation but uses the lookup structure
# that was built from the only initial pass
# over the GFF annotation file
# gene: Gene structure
# seg_lookup_table: Segment lookup table based on Gene ID
# strain_idx_table: Look-up strain index based on donor ID.
# segment_expr_info: Segment expression information
# gene_cds_begin_dict: Look-up to retrieve CDS beginnings associated with a gene [dict]
# ref_seq: Nucleotide sequence of chromosome associated with gene [str]
# mut_seq: Mutated sequence of donor
# fa_ptr: File handle to the output FASTA file
# mutation_mode: Mutation mode in {both, germline_only, somatic_only, None}
# size_factor: the adjusted weight from libsize
def annotate_gene_opt(gene=None, ref_seq=None, gene_idx=None,
                      seg_lookup_table=None, edge_lookup_table=None, size_factor=None, junction_list=None,
                      segment_expr_info=None, edge_expr_info=None, transcript_to_cds_table=None,
                      gene_to_transcript_table=None,
                      mutation_mode=None, mutation_sub_dic_vcf=None, mutation_sub_dic_maf=None, debug=False):


    sg = gene.splicegraph
    gene.from_sparse()

    output_id = 0
    output_peptide_list = []
    output_metadata_list = []

    # apply germline mutation
    # when germline mutation is applied, background_seq != ref_seq
    # otherwise, background_seq = ref_seq
    ref_mut_seq = apply_germline_mutation(ref_sequence=ref_seq, pos_start=gene.start, pos_end=gene.stop,
                                              mutation_sub_dic_vcf=mutation_sub_dic_vcf)

    # apply somatic mutation
    # som_exp_dict: (mutation_position) |-> (expression)
    # exon_som_dict: (exon_id) |-> (mutation_postion)
    som_exp_dict, exon_som_dict = None,None
    if mutation_sub_dic_maf is not None:
        exon_som_dict = get_exon_som_dict(gene, mutation_sub_dic_maf.keys())
        if segment_expr_info is not None:
            som_exp_dict = get_som_expr_dict(gene, mutation_sub_dic_maf.keys(), segment_expr_info, seg_lookup_table)

    # find background peptide
    # if no germline mutation is applies, germline key still exists, equals to reference.
    # return the list of the background peptide for each transcript
    background_pep_list = find_background_peptides(gene, ref_mut_seq['background'], gene_to_transcript_table,
                                                   transcript_to_cds_table)
    # check whether the junction (specific combination of vertices) also is annotated
    # as a junction of a protein coding transcript
    junction_flag = junction_is_annotated(gene, gene_to_transcript_table, transcript_to_cds_table)

    for idx in gene.vertex_order:
        n_read_frames = len(sg.reading_frames[idx])
        if n_read_frames == 0:  # no cds start, skip the vertex
            continue
        if is_isolated_cds(gene, idx):  # if it is an isolated cds, we add a flag idx '.', translate and output it
            gene.vertex_succ_list[idx].append('.')
        for prop_vertex in gene.vertex_succ_list[idx]:
            mut_seq_comb = get_mut_comb(exon_som_dict, idx, prop_vertex)
            for variant_comb in mut_seq_comb:  # go through each variant combination
                # Skip de-generate exons that contain less than one codon
                if gene.vertex_len_dict[idx] < 3:
                    continue
                for read_frame_tuple in sorted(sg.reading_frames[idx]):
                    if debug:
                        print(idx, prop_vertex, variant_comb, read_frame_tuple)
                    peptide_weight = 1.0 / n_read_frames
                    if prop_vertex != '.':
                        cross_peptide_mut, cross_peptide_ref, \
                        start_v1, stop_v1, start_v2, stop_v2, \
                        has_stop_codon, is_isolated, next_reading_frame = cross_peptide_result(read_frame_tuple, gene.strand, variant_comb, mutation_sub_dic_maf, ref_mut_seq, sg.vertices[:, prop_vertex])
                        if not has_stop_codon:
                            sg.reading_frames[prop_vertex].add(next_reading_frame)
                    else:
                        cross_peptide_mut, cross_peptide_ref, \
                        start_v1, stop_v1, start_v2, stop_v2, \
                        has_stop_codon, is_isolated = isolated_peptide_result(read_frame_tuple, gene.strand, variant_comb, mutation_sub_dic_maf,ref_mut_seq)

                    # If cross junction peptide has a stop-codon in it, the frame
                    # will not be propagated because the read is truncated before it reaches the end of the exon.
                    # also in mutation mode, only output the case where ref is different from mutated
                    if cross_peptide_mut != cross_peptide_ref or mutation_mode == 'ref':
                        match_ts_list = peptide_match(background_pep_list, cross_peptide_mut)
                        peptide_is_annotated = len(match_ts_list)
                        if not is_isolated:
                            junction_anno_flag = int(junction_flag[idx, prop_vertex])
                            if junction_list is not None:
                                if gene.strand == '+':
                                    junctionOI_flag = is_in_junction_list(sg.vertices[:, idx], sg.vertices[:, prop_vertex], gene.strand, junction_list)
                                else:
                                    junctionOI_flag = is_in_junction_list(sg.vertices[:, prop_vertex], sg.vertices[:, idx], gene.strand, junction_list)
                            else:
                                junctionOI_flag = '.'
                        else:
                            junction_anno_flag = '.'
                            junctionOI_flag = '.'
                        # Write the variant gene into the FASTA FP together with the donor ID
                        str_variant_comb = [str(ipos) for ipos in variant_comb]
                        if variant_comb != '.' and som_exp_dict is not None:  # which means there do exist some mutation
                            seg_exp_variant_comb = [str(som_exp_dict[ipos]) for ipos in variant_comb]
                        else:
                            seg_exp_variant_comb = '.'  # if no mutation, the segment expression is .
                        meta_header_line = "\t".join(
                            [str(gene_idx) + '.' + str(output_id), str(read_frame_tuple[2]), gene.name, gene.chr, gene.strand,
                             mutation_mode, "{:.3f}".format(peptide_weight), str(peptide_is_annotated), str(junction_anno_flag),
                             str(int(has_stop_codon)), str(junctionOI_flag), str(int(is_isolated)),
                             ';'.join(str_variant_comb), ';'.join(seg_exp_variant_comb)])
                        is_output = True
                        meta_header_line += ('\t' + str(start_v1) + ";" + str(stop_v1))
                        meta_header_line += (';' + str(start_v2) + ";" + str(stop_v2) + '\t')
                        meta_header_line += (str(idx) + ',' + str(prop_vertex) + '\t')

                        # deal with expression data
                        if edge_lookup_table is not None and not is_isolated:
                            sorted_pos = sp.sort(np.array([start_v1, stop_v1, start_v2,stop_v2]))
                            edge_expr = search_edge_metadata_segmentgraph(gene, sorted_pos, edge_lookup_table, edge_expr_info)
                            # edge_expr = edge_expr*size_factor
                        else:
                            edge_expr = '.'
                        meta_header_line += ("\t" .join([str(edge_expr)]))
                        output_metadata_list.append(meta_header_line)
                        peptide_str_pretty = '>' + str(gene_idx) + '.' + str(output_id) + '\n' + cross_peptide_mut
                        output_peptide_list.append(peptide_str_pretty)
                        output_id += 1
    if not sg.edges is None:
        gene.to_sparse()

    gene.processed = True
    return output_peptide_list,output_metadata_list
