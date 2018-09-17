from __future__ import print_function
# External libraries
import numpy as np
import scipy as sp
# immuno module
from immuno_filter import junction_is_annotated, peptide_match, find_background_peptides
from immuno_mutation import apply_germline_mutation,get_exon_som_dict,get_som_expr_dict,get_mut_comb
from utils import cross_peptide_result,is_isolated_cds,isolated_peptide_result,is_in_junction_list,get_segment_expr
from immuno_preprocess import search_edge_metadata_segmentgraph
from constant import NOT_EXIST

# Optimized annotation code that does not loop over the annotation but uses the lookup structure
# that was built from the only initial pass
# over the GFF annotation file
def calculate_output_peptide(gene=None, ref_seq=None, idx = None,
                      segments=None, edges=None, table=None,debug=False,size_factor=None, junction_list=None,
                      mutation = None):


    sg = gene.splicegraph
    gene.from_sparse()
    total_expr = 0

    output_id = 0
    output_peptide_list = []
    output_metadata_list = []

    # apply germline mutation
    # when germline mutation is applied, background_seq != ref_seq
    # otherwise, background_seq = ref_seq
    ref_mut_seq = apply_germline_mutation(ref_sequence=ref_seq, pos_start=gene.start, pos_end=gene.stop,
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
    background_pep_list = find_background_peptides(gene, ref_mut_seq['background'], table.gene_to_ts, table.ts_to_cds)

    # check whether the junction (specific combination of vertices) also is annotated
    # as a junction of a protein coding transcript
    junction_flag = junction_is_annotated(gene, table.gene_to_ts, table.ts_to_cds)
    reading_frame_dict = dict(sg.reading_frames)
    for v_id in gene.vertex_order:
        n_read_frames = len(reading_frame_dict[v_id])
        if n_read_frames == 0:  # no cds start, skip the vertex
            continue
        if is_isolated_cds(gene, v_id):  # if it is an isolated cds, we add a flag NOT_EXIST, translate and output it
            gene.vertex_succ_list[v_id].append(NOT_EXIST)
        for prop_vertex in gene.vertex_succ_list[v_id]:
            mut_seq_comb = get_mut_comb(exon_som_dict, v_id, prop_vertex)
            for variant_comb in mut_seq_comb:  # go through each variant combination
                # Skip de-generate exons that contain less than one codon
                if gene.vertex_len_dict[v_id] < 3:
                    continue
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
                    if peptide.mut != peptide.ref or mutation.mode == 'ref':
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
                            segment_expr = get_segment_expr(gene, coord, segments, idx)
                        else:
                            segment_expr = NOT_EXIST

                        meta_header_line = "\t".join(
                            [str(idx.gene) + '.' + str(output_id), str(read_frame_tuple[2]), gene.name, gene.chr, gene.strand,
                             mutation.mode, "{:.3f}".format(peptide_weight), str(peptide_is_annotated), str(junction_anno_flag),
                             str(int(flag.has_stop)), str(junctionOI_flag), str(int(flag.is_isolated)),
                             ';'.join(str_variant_comb), ';'.join(seg_exp_variant_comb)])

                        meta_header_line += ('\t' + str(coord.start_v1) + ";" + str(coord.stop_v1))
                        meta_header_line += (';' + str(coord.start_v2) + ";" + str(coord.stop_v2) + '\t')
                        meta_header_line += (str(v_id) + ',' + str(prop_vertex) + '\t')

                        # deal with expression data
                        if edges is not None and not flag.is_isolated:
                            sorted_pos = sp.sort(np.array([coord.start_v1, coord.stop_v1, coord.start_v2, coord.stop_v2]))
                            edge_expr = search_edge_metadata_segmentgraph(gene, sorted_pos, edges, idx)
                            total_expr += edge_expr
                            #edge_expr = edge_expr*size_factor
                        else:
                            edge_expr = NOT_EXIST
                        meta_header_line += ("\t" .join([str(edge_expr)]))
                        meta_header_line += "\t"+str(segment_expr)

                        output_metadata_list.append(meta_header_line)
                        peptide_str_pretty = '>' + str(idx.gene) + '.' + str(output_id) + '\n' + peptide.mut
                        output_peptide_list.append(peptide_str_pretty)
                        output_id += 1
    if not sg.edges is None:
        gene.to_sparse()

    gene.processed = True
    return output_peptide_list,output_metadata_list,total_expr
