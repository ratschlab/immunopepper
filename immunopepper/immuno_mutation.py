import bisect
from utils import get_all_comb


def apply_germline_mutation(ref_sequence, pos_start, pos_end, mutation_sub_dic_vcf):
    """
    Apply all the germline mutations to the reference sequence.


    Returns
    -------
    output_seq: dict, two keys. "ref" to the original reference sequence
                                "background" to the germline_mutation_applied ref (if there is germline mutation exist)
    """
    output_seq = {}
    output_seq['ref'] = ref_sequence  # copy the reference
    if mutation_sub_dic_vcf is not None:
        mut_seq = list(ref_sequence) # if no mutation infomation is provided
        variant_pos = [pos for pos in mutation_sub_dic_vcf.keys()]
        variant_pos_candi = [ipos for ipos in variant_pos if ipos > pos_start and ipos < pos_end]

        for variant_ipos in variant_pos_candi:
            ref_base = mutation_sub_dic_vcf[variant_ipos]['ref_base']
            mut_base = mutation_sub_dic_vcf[variant_ipos]['mut_base']
            mut_seq[variant_ipos] = mut_base
        output_seq['background'] = ''.join(mut_seq)
    else:
        output_seq['background'] = ref_sequence
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
        maf_path = args.maf_path[0]
        flag += 2
    map_dict = {0: 'ref', 1: 'germline', 2: 'somatic', 3: 'somatic_and_germline'}
    mutation_mode = map_dict[flag]
    return mutation_mode, vcf_path, maf_path


def get_exon_som_dict(gene, mutation_pos):
    """
    Build exon id to somatic mutation position dictionary map.
    """
    exon_list = gene.splicegraph.vertices
    exon_som_dict = {k:[] for k in range(exon_list.shape[1])}
    for ipos in mutation_pos:
        exon_id = bisect.bisect(exon_list[0, :], ipos)
        if exon_id > 0 and ipos <= exon_list[1][exon_id-1]:  # the mutation is within the pos
            exon_som_dict[exon_id-1].append(ipos)
    exon_som_dict['.'] = []  # for single cds case
    return exon_som_dict


def get_som_expr_dict(gene, mutation_pos, expr_mat, gene_ids_seg):
    """
    Build somatic mutation position to expression data dict map.
    """
    seg_mat = gene.segmentgraph.segments[0]
    som_expr_dict = {}
    seg_pos_list = gene_ids_seg[gene.name]
    for ipos in mutation_pos:
        seg_id = bisect.bisect(seg_mat,ipos)
        if seg_id > 0 and  ipos <= gene.segmentgraph.segments[1][seg_id-1]: # the mutation is within the pos
            expr = expr_mat[seg_pos_list[seg_id-1]]
            som_expr_dict[ipos] = expr
    return som_expr_dict


def get_mut_comb(exon_som_dict, idx, prop_vertex):
    """
    Get all the mutation combination given the mutation given.
    Parameters
    ----------
    exon_som_dict: dict, keys (exon id), values (somatic mutation position)
    idx: int, first exon id
    prop_vertex: int, second exon id

    Returns
    -------
    mut_comb: list of tuple, list of mutation combination

    """
    mut_comb = ['.']
    if exon_som_dict is not None:
        all_comb = get_all_comb(exon_som_dict[idx] + exon_som_dict[prop_vertex])
        mut_comb += all_comb
    return mut_comb

