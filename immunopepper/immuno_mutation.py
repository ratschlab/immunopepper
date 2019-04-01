"""Contain functions to deal with mutation"""
import bisect
import sys
from collections import namedtuple

from .constant import NOT_EXIST
from .immuno_preprocess import parse_mutation_from_maf,parse_mutation_from_vcf
from .utils import get_all_comb

import numpy as np

Mutation = namedtuple('Mutation', ['vcf_dict', 'maf_dict', 'mode'])


def apply_germline_mutation(ref_sequence, pos_start, pos_end, mutation_sub_dic_vcf):
    """Apply germline mutation on the reference sequence

    Parameters
    ----------
    ref_sequence: str. reference sequence of certain chromosome.
    pos_start: int. start position of sequence for applying germiline mutation.
    pos_end: int. Ending position of sequence for applying germline mutation.
    mutation_sub_dic_vcf: dict. (position) -> variant details

    Returns
    -------
    output_seq: dict. (sequence_type) -> list[char]. output['ref'] is the original
    reference sequence output['background'] is the germline-mutation-applied sequence
    if .maf file (somatic mutation) exists while is original sequence if no somatic
    information is available.
    """
    output_seq = {}
    output_seq['ref'] = ref_sequence  # copy the reference
    if mutation_sub_dic_vcf is not None:
        mut_seq = construct_mut_seq_with_str_concat(ref_sequence, pos_start, pos_end, mutation_sub_dic_vcf)
        output_seq['background'] = mut_seq
    else:
        output_seq['background'] = ref_sequence
    return output_seq


def apply_somatic_mutation(ref_sequence, pos_start, pos_end, mutation_sub_dic_maf):
    """Apply somatic mutation on the reference sequence

    Parameters
    ----------
    ref_sequence: str. reference sequence of certain chromosome.
    pos_start: int. start position of sequence for applying somatic mutation.
    pos_end: int. Ending position of sequence for applying somatic mutation.
    mutation_sub_dic_maf: dict. (position) -> variant details

    Returns
    -------
    output_seq: dict. (sequence_type) -> list[char]. output['ref'] is the original
    reference sequence output['background'] is the germline-mutation-applied sequence
    if .maf file (somatic mutation) exists while is original sequence if no somatic
    information is available.
    """
    if mutation_sub_dic_maf is not None:
        mut_seq = construct_mut_seq_with_str_concat(ref_sequence, pos_start, pos_end, mutation_sub_dic_maf)
    else:
        mut_seq = ref_sequence
    return mut_seq


def construct_mut_seq_with_str_concat(ref_seq, pos_start, pos_end, mut_dict):
    """ Applying germline mutation on given range and get mutated sequence

    Parameters
    ----------
    ref_seq: str. reference sequence
    pos_start: int. Start position of applying germline mutation
    pos_end: int. Stop position of applying germline mutation
    mut_dict: dict. Germline mutation dictioanry

    Returns
    -------
    mut_seq: str. mutation sequence
    """
    variant_pos_candi = [ipos for ipos in list(mut_dict.keys()) if ipos >= pos_start and ipos < pos_end]
    if len(variant_pos_candi) > 0:
        variant_pos_sorted = np.sort(variant_pos_candi)
        mut_seq_list = [ref_seq[:variant_pos_sorted[0]]]
        for i in range(len(variant_pos_sorted)-1):  # process all the mutation in order except the last one
            mut_base = mut_dict[variant_pos_sorted[i]]['mut_base']
            ref_base = mut_dict[variant_pos_sorted[i]]['ref_base']
            if mut_base != '*':
                mut_seq_list.append(mut_base)
            else:
                mut_seq_list.append(ref_base)
            mut_seq_list.append(ref_seq[variant_pos_sorted[i]+1:variant_pos_sorted[i+1]])
        # process the last mutation separately
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
    output_dir = args.output_dir
    is_error = True
    if mutation_mode == 'somatic_and_germline':
        if maf_file_path != '' and vcf_file_path != '':
            mutation_dic_maf = parse_mutation_from_maf(maf_file_path,output_dir)
            mutation_dic_vcf = parse_mutation_from_vcf(vcf_file_path,args.samples, args.heter_code)
            is_error = False
    elif mutation_mode == 'germline':
        if vcf_file_path != '':
            mutation_dic_maf = {}  # empty dic
            mutation_dic_vcf = parse_mutation_from_vcf(vcf_file_path,args.samples, args.heter_code)
            is_error = False
    elif mutation_mode == 'somatic':
        if maf_file_path != '':
            mutation_dic_maf = parse_mutation_from_maf(maf_file_path,output_dir)
            mutation_dic_vcf = {}
            is_error = False
    elif mutation_mode == 'ref':
        mutation_dic_maf = {}
        mutation_dic_vcf = {}
        is_error = False
    else:
        print('Mutation mode "%s" not recognized' % mutation_mode)

    if is_error:
         print("The input mutation file does not match the mutation mode (somatic, germline, somatic_and_germline), please check again")
         sys.exit(1)
    mutation = Mutation(mutation_mode,mutation_dic_maf,mutation_dic_vcf)
    return mutation


def get_exon_som_dict(gene, mutation_pos):
    """
    Build exon id (key) to somatic mutation position dictionary.
    """
    exon_list = gene.splicegraph.vertices
    exon_som_dict = {k:[] for k in range(exon_list.shape[1])}
    for ipos in mutation_pos:
        for i in range(exon_list.shape[1]):
            if ipos in range(exon_list[0,i], exon_list[1,i]):
                exon_som_dict[i].append(ipos)
    exon_som_dict[NOT_EXIST] = []  # for single cds case
    return exon_som_dict


def get_som_expr_dict(gene, mutation_pos, segments, Idx):
    """
    Build somatic mutation position(key) to expression data(value) dictionary.
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
    """ Get sub mutation namedtuple on given sample and chromosome """
    if (sample, chrm) in list(mutation.vcf_dict.keys()):
        mutation_sub_dict_vcf = mutation.vcf_dict[(sample, chrm)]
    else:
        mutation_sub_dict_vcf = None
    if (sample, chrm) in list(mutation.maf_dict.keys()):
        mutation_sub_dict_maf = mutation.maf_dict[(sample, chrm)]
    else:
        mutation_sub_dict_maf = None
    submutation = Mutation(mutation_sub_dict_vcf,mutation_sub_dict_maf,mutation.mode)
    return submutation
