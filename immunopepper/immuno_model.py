"""Contains all the output computation based on gene splicegraph"""


import numpy as np
import scipy as sp

from immunopepper.immuno_filter import junction_is_annotated, peptide_match,get_junction_anno_flag,get_filtered_metadata_list,get_full_peptide
from immunopepper.immuno_mutation import apply_germline_mutation,get_exon_som_dict,get_som_expr_dict,get_mut_comb,apply_somatic_mutation
from immunopepper.utils import cross_peptide_result,is_isolated_cds,isolated_peptide_result,is_in_junction_list,get_segment_expr,get_peptide_result,convert_namedtuple_to_str,write_namedtuple_list,get_total_gene_expr
from immunopepper.immuno_preprocess import search_edge_metadata_segmentgraph
from immunopepper.constant import NOT_EXIST
from immunopepper.immuno_nametuple import Output_metadata, Output_junc_peptide, Output_kmer, Simple_metadata, init_part_coord,Output_background

def get_simple_metadata(gene=None, ref_seq=None, idx=None,mutation=None, option=None):
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
       size_factor: Scalar. To adjust the expression counts based on the external file `libsize.tsv`
       junction_list: List. Work as a filter to indicate some exon pair has certain
           ordinary intron which can be ignored further.
       mutation: Namedtuple Mutation, store the mutation information of specific chromosome and sample.
           has the attribute ['mode', 'maf_dict', 'vcf_dict']

       Returns
       -------
       final_simple_meta:
       ref_mut_seq:
       exon_som_dict:

       """

    sg = gene.splicegraph
    gene.from_sparse()

    output_id = 0


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
        exon_som_dict = get_exon_som_dict(gene, list(mutation.maf_dict.keys()))

    simple_metadata_list = []
    reading_frame_dict = dict(sg.reading_frames)
    for v_id in gene.vertex_order:
        n_read_frames = len(reading_frame_dict[v_id])
        if n_read_frames == 0:  # no cds start, skip the vertex
            continue
        if len(gene.vertex_succ_list[v_id]) == 0: # if no successive vertex, we add a flag NOT_EXIST, translate and output it
            gene.vertex_succ_list[v_id].append(NOT_EXIST)
        for prop_vertex in gene.vertex_succ_list[v_id]:
            vertex_list = [v_id,prop_vertex]
            mut_seq_comb = get_mut_comb(exon_som_dict,vertex_list)
            for read_frame_tuple in sorted(reading_frame_dict[v_id]):
                has_stop_flag = True
                for variant_comb in mut_seq_comb:  # go through each variant combination
                    if option.debug:
                        print(v_id, prop_vertex, variant_comb, read_frame_tuple)
                    if prop_vertex != NOT_EXIST:
                        peptide, coord, flag, next_reading_frame = cross_peptide_result(read_frame_tuple, gene.strand, variant_comb, mutation.maf_dict, ref_mut_seq, sg.vertices[:, prop_vertex])
                        if not flag.has_stop:
                            reading_frame_dict[prop_vertex].add(next_reading_frame)
                    else:

                        peptide, coord, flag = isolated_peptide_result(read_frame_tuple, gene.strand, variant_comb, mutation.maf_dict,ref_mut_seq)
                    has_stop_flag = has_stop_flag and flag.has_stop
                gene_outputid = str(idx.gene) + '.' + str(output_id)
                simple_metadata = Simple_metadata(output_id=gene_outputid,
                                                  read_frame=read_frame_tuple,
                                                  exons_coor=coord,
                                                  vertex_idx=vertex_list,
                                                  has_stop_codon=has_stop_flag,
                                                  peptide_weight="{:.3f}".format(1/n_read_frames))
                simple_metadata_list.append(simple_metadata)
                output_id += 1

    concat_simple_meta_list = get_concat_metadata(gene, simple_metadata_list, option.kmer)
    if option.filter_redundant:
        filtered_simple_meta_list = get_filtered_metadata_list(simple_metadata_list,gene.strand)
    else:
        filtered_simple_meta_list = simple_metadata_list
    final_simple_meta = filtered_simple_meta_list+concat_simple_meta_list
    return final_simple_meta,ref_mut_seq,exon_som_dict

def get_and_write_peptide_and_kmer(gene=None, final_simple_meta=None, background_pep_list=None,ref_mut_seq=None, idx=None,
                         exon_som_dict=None,segments=None, edges=None, mutation=None,table=None,option=None,
                         size_factor=None, junction_list=None, filepointer=None):
    """

    Parameters
    ----------
    gene: Object, returned by SplAdder.
    final_simple_meta: Str, reference sequnce of specific chromosome
    idx: Namedtuple Idx, has attribute idx.gene and idx.sample
    segments: Namedtuple Segments, store segment expression information from count.hdf5.
       has attribute ['expr', 'lookup_table'].
    edges: Namedtuple Edges, store edges expression information from count.hdf5.
       has attribute ['expr','lookup_table']
    table: Namedtuple GeneTable, store the gene-transcript-cds mapping tables derived
       from .gtf file. has attribute ['gene_to_cds_begin', 'ts_to_cds', 'gene_to_cds']
    size_factor: Scalar. To adjust the expression counts based on the external file `libsize.tsv`
    junction_list: List. Work as a filter to indicate some exon pair has certain
       ordinary intron which can be ignored further.
    mutation: Namedtuple Mutation, store the mutation information of specific chromosome and sample.
       has the attribute ['mode', 'maf_dict', 'vcf_dict']
    filepointer:

    Returns
    -------

    """
    # check whether the junction (specific combination of vertices) also is annotated
    # as a  junction of a protein coding transcript
    junction_flag = junction_is_annotated(gene, table.gene_to_ts, table.ts_to_cds)
    total_expr = get_total_gene_expr(gene,segments,idx)
    som_exp_dict = get_som_expr_dict(gene, list(mutation.maf_dict.keys()), segments, idx)
    meta_field_list = ['output_id', 'read_frame', 'gene_name', 'gene_chr', 'gene_strand', 'mutation_mode',
                       'peptide_weight', 'peptide_annotated',
                       'junction_annotated', 'has_stop_codon', 'is_in_junction_list', 'is_isolated', 'variant_comb',
                       'variant_seg_expr',
                       'exons_coor', 'vertex_idx', 'junction_expr', 'segment_expr']
    junc_pep_field_list = ['output_id', 'id', 'new_line', 'peptide']
    kmer_field_list = ['kmer', 'id', 'expr', 'is_cross_junction']

    for simple_metadata in final_simple_meta:
        mut_seq_comb = get_mut_comb(exon_som_dict,simple_metadata.vertex_idx)
        gene_outputid = simple_metadata.output_id
        coord = simple_metadata.exons_coor
        vertex_list = simple_metadata.vertex_idx
        read_frame_tuple = simple_metadata.read_frame
        peptide_weight = simple_metadata.peptide_weight
        variant_id = 0
        for variant_comb in mut_seq_comb:  # go through each variant combination
            peptide,flag = get_peptide_result(simple_metadata, gene.strand, variant_comb, mutation.maf_dict, ref_mut_seq)

            # If cross junction peptide has a stop-codon in it, the frame
            # will not be propagated because the read is truncated before it reaches the end of the exon.
            # also in mutation mode, only output the case where ref is different from mutated
            if peptide.mut != peptide.ref or mutation.mode == 'ref' or option.output_silence:
                new_output_id = gene_outputid+'.'+str(variant_id)
                detail_id = gene.name+'.'+'_'.join([str(v) for v in vertex_list])+'.'+str(variant_id)
                match_ts_list = peptide_match(background_pep_list, peptide.mut)
                peptide_is_annotated = len(match_ts_list)
                junction_anno_flag = get_junction_anno_flag(junction_flag, vertex_list)
                junctionOI_flag = is_in_junction_list(gene.splicegraph,vertex_list, gene.strand, junction_list) # need re-write

                if variant_comb != NOT_EXIST and som_exp_dict is not None:  # which means there do exist some mutation
                    seg_exp_variant_comb = [int(som_exp_dict[ipos]) for ipos in variant_comb]
                else:
                    seg_exp_variant_comb = NOT_EXIST  # if no mutation or no count file,  the segment expression is .

                if segments is not None:
                    segment_expr, expr_list = get_segment_expr(gene, coord, segments, idx)
                else:
                    segment_expr, expr_list = NOT_EXIST, NOT_EXIST


                # deal with expression data
                if edges is not None and not flag.is_isolated:
                    edge_expr = search_edge_metadata_segmentgraph(gene, coord, edges, idx)
                    #edge_expr = edge_expr*size_factor
                else:
                    edge_expr = NOT_EXIST
                output_metadata = Output_metadata(output_id=new_output_id,
                                                  read_frame=read_frame_tuple.read_phase,
                                                  gene_name=gene.name,
                                                  gene_chr=gene.chr,
                                                  gene_strand=gene.strand,
                                                  mutation_mode=mutation.mode,
                                                  peptide_weight=peptide_weight,
                                                  peptide_annotated=peptide_is_annotated,
                                                  junction_annotated=junction_anno_flag,
                                                  has_stop_codon=int(flag.has_stop),
                                                  is_in_junction_list=junctionOI_flag,
                                                  is_isolated=int(flag.is_isolated),
                                                  variant_comb=variant_comb,
                                                  variant_seg_expr=seg_exp_variant_comb,
                                                  exons_coor=coord,
                                                  vertex_idx=vertex_list,
                                                  junction_expr=edge_expr,
                                                  segment_expr=segment_expr
                )
                variant_id += 1
                output_peptide = Output_junc_peptide(output_id='>'+new_output_id,
                                                id=detail_id,
                                                peptide=peptide.mut,exons_coor=coord)

                if option.kmer > 0:
                    output_kmer_list = create_output_kmer(output_peptide,expr_list,option.kmer)
                    write_namedtuple_list(filepointer.junction_kmer_fp, output_kmer_list, kmer_field_list)

                filepointer.junction_meta_fp.write(convert_namedtuple_to_str(output_metadata,meta_field_list)+'\n')
                filepointer.junction_peptide_fp.write(convert_namedtuple_to_str(output_peptide,junc_pep_field_list)+'\n')

    if not gene.splicegraph.edges is None:
        gene.to_sparse()
    gene.processed = True
    return total_expr

def get_and_write_background_peptide_and_kmer(gene, ref_mut_seq, table, Segments, Idx,filepointer,option):
    """Calculate the peptide translated from the complete transcript instead of single exon pairs

    Parameters
    ----------
    gene: Object. Created by SplAdder
    ref_seq: List(str). Reference sequence of certain chromosome.
    table: Namedtuple GeneTable, store the gene-transcript-cds mapping tables derived
       from .gtf file. has attribute ['gene_to_cds_begin', 'ts_to_cds', 'gene_to_cds']

    Returns
    -------
    peptide_list: List[str]. List of all the peptide translated from the given
       splicegraph and annotation.
    (ts_list): List[str]. List of all the transcript indicated by the  annotation file
        can be used to generate artifical reads.
    """
    ref_seq = ref_mut_seq['background']
    gene_to_transcript_table,transcript_cds_table = table.gene_to_ts, table.ts_to_cds,
    gene_transcripts = gene_to_transcript_table[gene.name]
    back_pep_field_list = ['id', 'new_line', 'peptide']
    kmer_field_list = ['kmer', 'id', 'expr', 'is_cross_junction']
    background_peptide_list = []
    # Generate a background peptide for every variant transcript
    for ts in gene_transcripts:
        # No CDS entries for transcript in annotation file...
        if ts not in transcript_cds_table:
            #print("WARNING: Transcript not in CDS table")
            continue
        cds_list = transcript_cds_table[ts]
        cds_expr_list, cds_string, cds_peptide = get_full_peptide(gene,ref_seq,cds_list,Segments,Idx,mode='back')
        peptide = Output_background(ts,cds_peptide)
        background_peptide_list.append(peptide)
        filepointer.background_peptide_fp.write(convert_namedtuple_to_str(peptide,back_pep_field_list)+'\n')
        if option.kmer > 0:
            output_kmer_list = create_output_kmer(peptide, cds_expr_list, option.kmer)
            write_namedtuple_list(filepointer.background_kmer_fp,output_kmer_list,kmer_field_list)
    gene.processed = True
    return background_peptide_list


def create_output_kmer(output_peptide, expr_list, k):
    """Calculate the output kmer and the corresponding expression based on output peptide

    Parameters
    ----------
    peptide_list: Output_junc_peptide. Filtered output_peptide_list.
    expr_lists: List(Tuple). Filtered expr_list.
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
        m = int(L1 / 3)
        if L1%3 == 0:
            spanning_id_range = list(range(max(m-k+1,0),m))
        else:
            spanning_id_range = list(range(max(m-k+1,0),m+1))
        return spanning_id_range

    output_kmer_list = []
    peptide = output_peptide.peptide
    peptide_head = output_peptide.id
    if hasattr(output_peptide,'exons_coor'):
        coord = output_peptide.exons_coor
        spanning_index = get_spanning_index(coord,k)
    else:
        spanning_index = NOT_EXIST
    # decide the kmer that spans over the cross junction
    expr_array = change_expr_lists_to_array(expr_list)
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
            output_kmer_list.append(kmer)
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
        output_kmer_list.append(kmer)
    return output_kmer_list

def get_concat_metadata(gene,output_metadata_list,k):
    def in_the_same_read_frame(front_coord, back_coord, strand):
        if strand == '+':
            return (front_coord.stop_v2-back_coord.start_v1)%3 == 0
        else:
            return (back_coord.stop_v1-front_coord.start_v2)%3 == 0
    concat_simple_meta_list = []
    vertices = gene.splicegraph.vertices
    vertex_len = vertices[1,:]-vertices[0,:]
    key_id_list = np.where(vertex_len < (k+1)*3)[0]
    vertex_id_pair_list = [metadata.vertex_idx for metadata in output_metadata_list]
    stop_codon_list = [metadata.has_stop_codon for metadata in output_metadata_list]
    strand = gene.strand
    for key_id in key_id_list:
        front_id_list =[i for i, vert_pair in enumerate(vertex_id_pair_list) if vert_pair[1] == key_id
                        and not stop_codon_list[i]]
        back_id_list =[i for i, vert_pair in enumerate(vertex_id_pair_list) if vert_pair[0] == key_id
                       and vert_pair[1] != NOT_EXIST]
        for front_id in front_id_list:
            for back_id in back_id_list:
                front_meta = output_metadata_list[front_id]
                back_meta = output_metadata_list[back_id]
                front_coord = front_meta.exons_coor
                back_coord = back_meta.exons_coor
                if not in_the_same_read_frame(front_coord,back_coord,strand):
                    continue
                front_vertex_id = front_meta.vertex_idx
                middle_exon_coord = gene.splicegraph.vertices[:,front_vertex_id[1]]
                triple_coord = init_part_coord(start_v1=front_meta.exons_coor.start_v1,
                                               stop_v1=front_meta.exons_coor.stop_v1,
                                               start_v2=middle_exon_coord[0],
                                               stop_v2=middle_exon_coord[1],
                                               start_v3=back_meta.exons_coor.start_v2,
                                               stop_v3=back_meta.exons_coor.stop_v2)
                triple_output_id = front_meta.output_id+'_'+back_meta.output_id.split('.')[-1]
                triple_frame = front_meta.read_frame
                triple_vertex_idx = front_meta.vertex_idx+[back_meta.vertex_idx[-1]]
                triple_has_stop = back_meta.has_stop_codon
                new_simple_metadata = Simple_metadata(output_id=triple_output_id,
                                                      read_frame=triple_frame,
                                                      exons_coor=triple_coord,
                                                      vertex_idx=triple_vertex_idx,
                                                      has_stop_codon=triple_has_stop,
                                                      peptide_weight=front_meta.peptide_weight)
                concat_simple_meta_list.append(new_simple_metadata)
    return concat_simple_meta_list
