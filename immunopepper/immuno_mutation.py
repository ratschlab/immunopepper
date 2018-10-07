import bisect
import sys
from collections import namedtuple

from constant import NOT_EXIST
from immuno_preprocess import parse_mutation_from_maf,parse_mutation_from_vcf
from utils import mut_replace, get_all_comb


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
    output_seq['ref'] = list(ref_sequence)  # copy the reference
    if mutation_sub_dic_vcf is not None:
        mut_seq = list(ref_sequence) # if no mutation infomation is provided
        variant_pos = [pos for pos in mutation_sub_dic_vcf.keys()]
        variant_pos_candi = [ipos for ipos in variant_pos if ipos > pos_start and ipos < pos_end]

        for variant_ipos in variant_pos_candi:
            ref_base = mutation_sub_dic_vcf[variant_ipos]['ref_base']
            mut_base = mutation_sub_dic_vcf[variant_ipos]['mut_base']
            mut_seq = mut_replace(mut_seq, variant_ipos, ref_base,mut_base)
        output_seq['background'] = mut_seq
    else:
        output_seq['background'] = list(ref_sequence)
    return output_seq


def get_mutation_mode_from_parser(args):
    """

    Parameters
    ----------
    args

    Returns
    -------

    """
    Mutation = namedtuple('Mutation', ['mode','maf_dict','vcf_dict'])
    mutation_mode = args.mutation_mode
    maf_file_path = args.maf_path[0]
    vcf_file_path = args.vcf_path[0]
    is_error = True
    if mutation_mode == 'somatic_and_germline':
        if maf_file_path != '' and vcf_file_path != '':
            mutation_dic_maf = parse_mutation_from_maf(maf_file_path)
            mutation_dic_vcf = parse_mutation_from_vcf(vcf_file_path)
            is_error = False
    elif mutation_mode == 'germline':
        if vcf_file_path != '':
            mutation_dic_maf = {}  # empty dic
            mutation_dic_vcf = parse_mutation_from_vcf(vcf_file_path)
            is_error = False
    elif mutation_mode == 'somatic':
        if maf_file_path != '':
            mutation_dic_maf = parse_mutation_from_maf(maf_file_path)
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
