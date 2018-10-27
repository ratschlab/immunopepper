import bisect
from utils import get_all_comb
from collections import namedtuple
from constant import NOT_EXIST
from immuno_preprocess import parse_mutation_from_maf,parse_mutation_from_vcf
import sys

import numpy as np

Mutation = namedtuple('Mutation', ['vcf_dict', 'maf_dict', 'mode'])

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
        mut_seq = construct_mut_seq_with_str_concat(ref_sequence, pos_start, pos_end, mutation_sub_dic_vcf)
        output_seq['background'] = mut_seq
    else:
        output_seq['background'] = ref_sequence
    return output_seq



def construct_mut_seq_with_str_concat(ref_seq, pos_start, pos_end, mut_dict):
    variant_pos_candi = [ipos for ipos in mut_dict.keys() if ipos > pos_start and ipos < pos_end]
    if len(variant_pos_candi) >0 :
        variant_pos_sorted = np.sort(variant_pos_candi)
        mut_seq_list = [ref_seq[:variant_pos_sorted[0]]]
        for i in range(len(variant_pos_sorted)-1):
            mut_base = mut_dict[variant_pos_sorted[i]]['mut_base']
            ref_base = mut_dict[variant_pos_sorted[i]]['ref_base']
            if mut_base != '*':
                mut_seq_list.append(mut_base)
            else:
                mut_seq_list.append(ref_base)
            mut_seq_list.append(ref_seq[variant_pos_sorted[i]+1:variant_pos_sorted[i+1]])
        mut_base = mut_dict[variant_pos_sorted[-1]]['mut_base']
        ref_base = mut_dict[variant_pos_sorted[-1]]['ref_base']
        if mut_base != '*':
            mut_seq_list.append(mut_base)
        else:
            mut_seq_list.append(ref_base)
        mut_seq_list.append(ref_seq[variant_pos_sorted[i]+1:variant_pos_sorted[i+1]])
        mut_base = mut_dict[variant_pos_sorted[-1]]['mut_base']
        ref_base = mut_dict[variant_pos_sorted[-1]]['ref_base']
        if mut_base != '*':
            mut_seq_list.append(mut_base)
        else:
            mut_seq_list.append(ref_base)
        mut_seq_list.append(ref_seq[variant_pos_sorted[-1]+1:])
        mut_seq = ''.join(mut_seq_list)
    else:
        mut_seq = ref_seq
    return mut_seq

def get_mutation_mode_from_parser(args):
    Mutation = namedtuple('Mutation', ['mode','maf_dict','vcf_dict'])
    mutation_mode = args.mutation_mode
    maf_file_path = args.maf_path
    vcf_file_path = args.vcf_path
    is_error = True
    if mutation_mode == 'somatic_and_germline':
        if maf_file_path[0] != '' and vcf_file_path[0] != '':
            mutation_dic_maf = parse_mutation_from_maf(maf_file_path[0])
            mutation_dic_vcf = parse_mutation_from_vcf(vcf_file_path[0])
            is_error = False
    elif mutation_mode == 'germline':
        if vcf_file_path[0] != '':
            mutation_dic_maf = {}  # empty dic
            mutation_dic_vcf = parse_mutation_from_vcf(vcf_file_path[0])
            is_error = False
    elif mutation_mode == 'somatic':
        if maf_file_path[0] != '':
            mutation_dic_maf = parse_mutation_from_maf(maf_file_path[0])
            mutation_dic_vcf = {}
            is_error = False
    elif mutation_mode == 'ref':
        mutation_dic_maf = {}
        mutation_dic_vcf = {}
        is_error = False
    else:
        print('Mutation mode "%s" not recognized' % mutation_mode)

    if is_error:
         print("The input mutation file does not match the mutation mode, please check again")
         sys.exit(1)
    mutation = Mutation(mutation_mode,mutation_dic_maf,mutation_dic_vcf)
    return mutation


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
    exon_som_dict[NOT_EXIST] = []  # for single cds case
    return exon_som_dict


def get_som_expr_dict(gene, mutation_pos, segments, Idx):
    """
    Build somatic mutation position to expression data dict map.
    """
    seg_mat = gene.segmentgraph.segments[0]
    som_expr_dict = {}

    seg_pos_list = segments.lookup_table[gene.name]
    for ipos in mutation_pos:
        seg_id = bisect.bisect(seg_mat,ipos)
        if seg_id > 0 and ipos <= gene.segmentgraph.segments[1][seg_id-1]: # the mutation is within the pos
            expr = segments.expr[seg_pos_list[seg_id-1],Idx.sample]
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
    mut_comb = [NOT_EXIST]
    if exon_som_dict is not None:
        all_comb = get_all_comb(exon_som_dict[idx] + exon_som_dict[prop_vertex])
        mut_comb += all_comb
    return mut_comb


def get_sub_mutation_tuple(mutation, sample, chrm):
    if (sample, chrm) in mutation.vcf_dict.keys():
        mutation_sub_dict_vcf = mutation.vcf_dict[(sample, chrm)]
    else:
        mutation_sub_dict_vcf = None
    if (sample, chrm) in mutation.maf_dict.keys():
        mutation_sub_dict_maf = mutation.maf_dict[(sample, chrm)]
    else:
        mutation_sub_dict_maf = None
    submutation = Mutation(mutation_sub_dict_vcf,mutation_sub_dict_maf,mutation.mode)
    return submutation
