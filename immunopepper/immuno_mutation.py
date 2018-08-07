import numpy as np
import bisect
from utils import get_all_comb

# if is_code_somatic = False, the maf information would not be encoded,
# even the mutation mode includes 'somatic'
def apply_germline_mutation(ref_sequence, pos_start,pos_end,mutation_sub_dic_vcf):
    ## Todo: actually we can apply slicing on the reference geen so that save the memory consumed
    ## Todo: split into two functions: germline and somatic
    output_seq = {}
    output_seq['ref'] = ref_sequence  # copy the reference
    mut_seq = list(ref_sequence) # if no mutation infomation is provided

    variant_pos = [pos for pos in mutation_sub_dic_vcf.keys()]
    variant_pos_candi = [ipos for ipos in variant_pos if ipos > pos_start and ipos < pos_end]

    for variant_ipos in variant_pos_candi:
        ref_base = mutation_sub_dic_vcf[variant_ipos]['ref_base']
        mut_base = mutation_sub_dic_vcf[variant_ipos]['mut_base']
       # print ref_sequence[variant_ipos]
       # print ref_base
        #assert ref_sequence[variant_ipos] == ref_base
        mut_seq[variant_ipos] = mut_base
    output_seq['germline'] = ''.join(mut_seq)
    return output_seq

def apply_somatic_mutation(ref_sequence,all_comb, mutation_sub_dic_maf):
    output_seq = {}
    output_seq['.'] = ref_sequence
    for icomb in all_comb:
        mut_seq = list(ref_sequence)
        # print(pos_start, pos_end, icomb)
        for variant_ipos in icomb:
            ref_base = mutation_sub_dic_maf[variant_ipos]['ref_base']
            mut_base = mutation_sub_dic_maf[variant_ipos]['mut_base']
            strand = mutation_sub_dic_maf[variant_ipos]['strand']
            variant_Classification = mutation_sub_dic_maf[variant_ipos]['variant_Classification']
            variant_Type = mutation_sub_dic_maf[variant_ipos]['variant_Type']
            # assert ref_sequence[variant_ipos] == ref_base
            mut_seq[variant_ipos] = mut_base
            output_seq[icomb] = ''.join(mut_seq)
    return output_seq

def get_mutation_mode_from_parser(args):
    flag = 0
    if len(args.vcf_path) == 0:
        vcf_path = ''
    else:
        vcf_path = args.vcf_path[0]
        flag += 1

    if len(args.maf_path) == 0:
        maf_path = ''
    else:
        maf_path = args.maf_path
        flag += 2
    map_dict = {0:'ref', 1:'germline', 2:'somatic', 3:'somatic_and_germline'}
    mutation_mode = map_dict[flag]
    return mutation_mode, vcf_path, maf_path

def get_exon_som_dict(gene,mutation_pos):
    seg_mat = gene.segmentgraph.segments[0]
    seg_exon_mat = gene.segmentgraph.seg_match
    exon_som_dict = {k:[] for k in range(seg_exon_mat.shape[0])}
    for ipos in mutation_pos:
        seg_id = bisect.bisect(seg_mat,ipos)
        if seg_id > 0 and  ipos <= gene.segmentgraph.segments[1][seg_id-1]: # the mutation is within the pos
            exon_list = np.where(seg_exon_mat[:,seg_id-1] == 1)[0]
            for id_exon in exon_list:
                exon_som_dict[id_exon].append(ipos)
    exon_som_dict['.'] = [] # for single cds case
    return exon_som_dict

def get_som_expr_dict(gene,mutation_pos,expr_mat,gene_ids_seg):
    seg_mat = gene.segmentgraph.segments[0]
    som_expr_dict = {}
    seg_pos_list = gene_ids_seg[gene.name]
    for ipos in mutation_pos:
        seg_id = bisect.bisect(seg_mat,ipos)
        if seg_id > 0 and  ipos <= gene.segmentgraph.segments[1][seg_id-1]: # the mutation is within the pos
            expr = expr_mat[seg_pos_list[seg_id-1]]
            som_expr_dict[ipos] = expr
    return som_expr_dict

def get_mut_seq_dict(background_seq, mutation_sub_dic_maf, exon_som_dict, idx, prop_vertex):
    """
    Get sequence list when applied all the mutation combinations

    Parameters
    ----------
    background_seq: str, background_seq
    mutation_sub_dic_maf: somatic mutation dict, {mutation_position:mutation_detail}
    exon_som_dict: dict, {exon_id:mutation_position}
    idx: int, vertex index
    prop_vertex: int, next vertex index

    Returns
    -------

    """
    mut_seq_dict = {'.': background_seq}
    if mutation_sub_dic_maf is not None:
        all_comb = get_all_comb(exon_som_dict[idx] + exon_som_dict[prop_vertex])
        if len(all_comb) > 0:
            mut_seq_dict = apply_somatic_mutation(ref_sequence=background_seq, all_comb=all_comb,
                                                  mutation_sub_dic_maf=mutation_sub_dic_maf)
    return mut_seq_dict


def get_mut_comb(mutation_sub_dic_maf, exon_som_dict, idx, prop_vertex):
    mut_comb = ['.']
    if mutation_sub_dic_maf is not None:
        all_comb = get_all_comb(exon_som_dict[idx] + exon_som_dict[prop_vertex])
        mut_comb += all_comb
    return mut_comb

