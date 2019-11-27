"""Contain functions to deal with mutation"""
import bisect
import sys
from collections import namedtuple
from functools import reduce
import logging

from .constant import NOT_EXIST
from .immuno_preprocess import parse_mutation_from_maf,parse_mutation_from_vcf
from .utils import get_all_comb

import numpy as np

Mutation = namedtuple('Mutation', ['mode','germline_mutation_dict','somatic_mutation_dict'])


def apply_germline_mutation(ref_sequence, pos_start, pos_end, mutation_sub_dict):
    """Apply germline mutation on the reference sequence

    Parameters
    ----------
    ref_sequence: str. reference sequence of certain chromosome.
    pos_start: int. start position of sequence for applying germiline mutation.
    pos_end: int. Ending position of sequence for applying germline mutation.
    mutation_sub_dict_vcf: dict. (position) -> variant details

    Returns
    -------
    output_seq: dict. (sequence_type) -> list[char]. output['ref'] is the original
    reference sequence output['background'] is the germline-mutation-applied sequence
    if .maf file (somatic mutation) exists while is original sequence if no somatic
    information is available.
    """
    output_seq = {}
    output_seq['ref'] = ref_sequence  # copy the reference
    if mutation_sub_dict is not None:
        mut_seq = construct_mut_seq_with_str_concat(ref_sequence, pos_start, pos_end, mutation_sub_dict)
        output_seq['background'] = mut_seq
    else:
        output_seq['background'] = ref_sequence
    return output_seq


def apply_somatic_mutation(ref_sequence, pos_start, pos_end, mutation_sub_dict):
    """Apply somatic mutation on the reference sequence

    Parameters
    ----------
    ref_sequence: str. reference sequence of certain chromosome.
    pos_start: int. start position of sequence for applying somatic mutation.
    pos_end: int. Ending position of sequence for applying somatic mutation.
    mutation_sub_dict: dict. (position) -> variant details

    Returns
    -------
    output_seq: dict. (sequence_type) -> list[char]. output['ref'] is the original
    reference sequence output['background'] is the germline-mutation-applied sequence
    if .maf file (somatic mutation) exists while is original sequence if no somatic
    information is available.
    """
    if mutation_sub_dict is not None:
        mut_seq = construct_mut_seq_with_str_concat(ref_sequence, pos_start, pos_end, mutation_sub_dict)
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

def parse_mutation_file(mutation_file_path,output_dir,heter_code,mut_pickle=False,h5_sample_list=None):
    if mutation_file_path.lower().endswith('.maf'):
        mutation_dict = parse_mutation_from_maf(maf_path=mutation_file_path,output_dir=output_dir,mut_pickle=mut_pickle)
    elif mutation_file_path.lower().endswith('.vcf') or mutation_file_path.lower().endswith('.h5'): # we also accept hdf5 file format
        mutation_dict = parse_mutation_from_vcf(vcf_path=mutation_file_path,output_dir=output_dir,mut_pickle=mut_pickle,
                                                heter_code=heter_code,h5_sample_list=h5_sample_list)
    else:
        logging.error("Invalid mutation files. Please ensure it is in maf or vcf format.")
        sys.exit(1)
    return mutation_dict


def get_mutation_mode_from_parser(args):
    "Check if the input files match the mutation mode"
    mutation_mode = args.mutation_mode
    germline_file_path = args.germline
    somatic_file_path = args.somatic
    output_dir = args.output_dir
    heter_code = args.heter_code
    mut_pickle = args.use_mut_pickle
    h5_sample_list = args.samples
    is_error = True
    if mutation_mode == 'somatic_and_germline':
        if somatic_file_path != '' and germline_file_path != '':
            somatic_mutation_dict = parse_mutation_file(somatic_file_path,output_dir,heter_code,mut_pickle,h5_sample_list)
            germline_mutation_dict = parse_mutation_file(germline_file_path,output_dir,heter_code,mut_pickle,h5_sample_list)
            is_error = False
    elif mutation_mode == 'germline':
        if germline_file_path != '':
            somatic_mutation_dict = {}  # empty dic
            germline_mutation_dict = parse_mutation_file(germline_file_path,output_dir,heter_code,mut_pickle,h5_sample_list)
            is_error = False
    elif mutation_mode == 'somatic':
        if somatic_file_path != '':
            somatic_mutation_dict = parse_mutation_file(somatic_file_path,output_dir,heter_code,mut_pickle,h5_sample_list)
            germline_mutation_dict = {}
            is_error = False
    elif mutation_mode == 'ref':
        somatic_mutation_dict = {}
        germline_mutation_dict = {}
        is_error = False
    else:
        logging.error('Mutation mode "%s" not recognized, please check again.' % mutation_mode)

    if is_error:
         logging.error("immuno_mutation.py: The input mutation file does not match the mutation mode (somatic, germline, somatic_and_germline), please check again")
         sys.exit(1)
    mutation = Mutation(mutation_mode,germline_mutation_dict=germline_mutation_dict,somatic_mutation_dict=somatic_mutation_dict)
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
    if segments is None:
        return None
    seg_mat = gene.segmentgraph.segments[0]
    som_expr_dict = {}

    seg_pos_list = segments.lookup_table[gene.name]
    for ipos in mutation_pos:
        seg_id = bisect.bisect(seg_mat,ipos)
        if seg_id > 0 and ipos <= gene.segmentgraph.segments[1][seg_id-1]: # the mutation is within the pos
            expr = segments.expr[seg_pos_list[seg_id-1],Idx.sample]
            som_expr_dict[ipos] = expr
    return som_expr_dict


def get_mut_comb(exon_som_dict,vertex_list):
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
        exon_list = map(lambda x: exon_som_dict[x],vertex_list)
        all_comb = get_all_comb(reduce(lambda x,y:x+y,exon_list))
        mut_comb += all_comb
    return mut_comb


def get_sub_mutation_tuple(mutation, sample, chrm):
    """ Get sub mutation namedtuple on given sample and chromosome """
    mode_map = {'ref':0,'somatic':1,'germline':2,'somatic_and_germline':3}
    inv_mode_map = {0:'ref',1:'somatic',2:'germline',3:'somatic_and_germline'}
    loaded_mode_code = 0
    expected_mode_code = mode_map[mutation.mode]
    if (sample, chrm) in list(mutation.germline_mutation_dict.keys()):
        germline_mutation_sub_dict = mutation.germline_mutation_dict[(sample, chrm)]
        loaded_mode_code += mode_map['germline']
    else:
        germline_mutation_sub_dict = {}
    if (sample, chrm) in list(mutation.somatic_mutation_dict.keys()):
        somatic_mutation_sub_dict = mutation.somatic_mutation_dict[(sample, chrm)]
        loaded_mode_code += mode_map['somatic']
    else:
        somatic_mutation_sub_dict = {}
    if expected_mode_code != loaded_mode_code:
        logging.warning("The expected mode is {} but the loaded mode is {}."
              " Probably there is no mutation in the given chromosome for the given sample.".format(inv_mode_map[expected_mode_code],inv_mode_map[loaded_mode_code]))
    submutation = Mutation(mode=mutation.mode,somatic_mutation_dict=somatic_mutation_sub_dict,germline_mutation_dict=germline_mutation_sub_dict)
    return submutation
