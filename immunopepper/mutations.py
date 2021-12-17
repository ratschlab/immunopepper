"""Contain functions to deal with mutation"""
import bisect
import logging
import sys
import numpy as np
import pysam
import h5py

from functools import reduce

from .namedtuples import Mutation
from .preprocess import parse_mutation_from_maf
from .preprocess import parse_mutation_from_vcf
from .utils import get_all_comb


def apply_germline_mutation(ref_sequence_file, chrm, pos_start, pos_end, mutation_sub_dict):
    """Apply germline mutation on the reference sequence

    Parameters
    ----------
    ref_sequence: str. reference sequence of certain chromosome.
    pos_start: int. start position of sequence for applying germiline mutation.
    pos_end: int. Ending position of sequence for applying germline mutation.
    mutation_sub_dict: dict. (position) -> variant details

    Returns
    -------
    output_seq: dict. (sequence_type) -> list[char]. where
        output['ref'] is the original reference sequence 
        output['background'] is the reference seq with germline-mutation-applied
    """
    with pysam.FastaFile(ref_sequence_file) as fh:
        ref_sequence = fh.fetch(chrm, pos_start, pos_end)
    output_seq = {}
    output_seq['ref'] = ref_sequence  # copy the reference
    if mutation_sub_dict is not None:
        mut_seq = construct_mut_seq_with_str_concat(ref_sequence, pos_start, pos_end, mutation_sub_dict)
        output_seq['background'] = mut_seq
    else:
        output_seq['background'] = ref_sequence
    return output_seq


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
        mut_seq_list = [ref_seq[:variant_pos_sorted[0] - pos_start]]
        for i in range(len(variant_pos_sorted)-1):  # process all the mutation in order except the last one
            mut_base = mut_dict[variant_pos_sorted[i]]['mut_base']
            ref_base = mut_dict[variant_pos_sorted[i]]['ref_base']
            if mut_base != '*':
                mut_seq_list.append(mut_base)
            else:
                mut_seq_list.append(ref_base)
            mut_seq_list.append(ref_seq[variant_pos_sorted[i] - pos_start + 1:variant_pos_sorted[i + 1] - pos_start])
        # process the last mutation separately
        mut_base = mut_dict[variant_pos_sorted[-1]]['mut_base']
        ref_base = mut_dict[variant_pos_sorted[-1]]['ref_base']
        if mut_base != '*':
            mut_seq_list.append(mut_base)
        else:
            mut_seq_list.append(ref_base)
        mut_seq_list.append(ref_seq[variant_pos_sorted[-1] -pos_start + 1 :])
        mut_seq = ''.join(mut_seq_list)
    else:
        mut_seq = ref_seq
    return mut_seq


def parse_mutation_file(mutation_tag, mutation_file_path, output_dir, heter_code, mut_pickle=False, target_sample_list=None, mutation_sample=None, name_eq_dict={}):
    if mutation_file_path.lower().endswith('.maf'):
        mutation_dict = parse_mutation_from_maf(mutation_tag=mutation_tag, target_sample_list=target_sample_list, mutation_sample=mutation_sample, maf_path=mutation_file_path,output_dir=output_dir,mut_pickle=mut_pickle,sample_eq_dict=name_eq_dict)
    elif mutation_file_path.lower().endswith('.vcf') or mutation_file_path.lower().endswith('.h5'): # we also accept hdf5 file format
        mutation_dict = parse_mutation_from_vcf(mutation_tag=mutation_tag, vcf_path=mutation_file_path,output_dir=output_dir,mut_pickle=mut_pickle,
                                                heter_code=heter_code,target_sample_list=target_sample_list, mutation_sample=mutation_sample, sample_eq_dict=name_eq_dict)
    else:
        logging.error("Invalid mutation files. Please ensure it is in maf or vcf format.")
        sys.exit(1)
    return mutation_dict


def get_mutation_mode_from_parser(args, target_sample_list):
    "Check if the input files match the mutation mode"
    mutation_mode = args.mutation_mode
    germline_file_path = args.germline
    somatic_file_path = args.somatic
    output_dir = args.output_dir
    heter_code = args.heter_code
    mut_pickle = args.use_mut_pickle
    is_error = True
    graph_to_somatic_names = {}
    graph_to_germline_names = {}
    sample_name_map_tbl = args.sample_name_map
    if sample_name_map_tbl is not None:
        with open(sample_name_map_tbl, "r") as f:
            for line in f.readlines():
                equivalence = line.split("\n")[0].split('\t')
                if len(equivalence) == 2:
                    graph_to_somatic_names[equivalence[0]] = equivalence[1]
                    graph_to_germline_names[equivalence[0]] = equivalence[1]
                elif len(equivalence) == 3:
                    graph_to_germline_names[equivalence[0]] = equivalence[1]
                    graph_to_somatic_names[equivalence[0]] = equivalence[2]
                else:
                    logging.error("--args.sample_name_map provided with wrong input format, please check documentation")
                    sys.exit(1)
    else:
        for sample in target_sample_list:
            graph_to_somatic_names[sample] = sample
        graph_to_germline_names = graph_to_somatic_names
    if mutation_mode == 'somatic_and_germline':
        if somatic_file_path != '' and germline_file_path != '':
            somatic_mutation_dict = parse_mutation_file("somatic", somatic_file_path, output_dir, heter_code, mut_pickle, target_sample_list, args.mutation_sample, graph_to_somatic_names)
            germline_mutation_dict = parse_mutation_file("germline", germline_file_path, output_dir, heter_code, mut_pickle, target_sample_list, args.mutation_sample, graph_to_germline_names)
            is_error = False
    elif mutation_mode == 'germline':
        if germline_file_path != '':
            somatic_mutation_dict = {}  # empty dic
            germline_mutation_dict = parse_mutation_file("germline", germline_file_path, output_dir, heter_code, mut_pickle, target_sample_list, args.mutation_sample, graph_to_germline_names)
            is_error = False
    elif mutation_mode == 'somatic':
        if somatic_file_path != '':
            somatic_mutation_dict = parse_mutation_file("somatic", somatic_file_path, output_dir, heter_code, mut_pickle, target_sample_list, args.mutation_sample, graph_to_somatic_names)
            germline_mutation_dict = {}
            is_error = False
    elif mutation_mode == 'ref':
        somatic_mutation_dict = {}
        germline_mutation_dict = {}
        is_error = False
    else:
        logging.error('Mutation mode "%s" not recognized, please check again.' % mutation_mode)

    if is_error:
         logging.error("mutations.py: The input mutation file does not match the mutation mode (somatic, germline, somatic_and_germline), please check again")
         sys.exit(1)
    mutation = Mutation(mutation_mode, germline_mutation_dict=germline_mutation_dict, somatic_mutation_dict=somatic_mutation_dict)
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
    exon_som_dict[np.nan] = []  # for single cds case
    return exon_som_dict


def get_som_expr_dict(gene, mutation_pos, countinfo, seg_counts, mut_count_id):
    """
    Build somatic mutation position(key) to expression data(value) dictionary.
    """
    if countinfo is None:
        return None
    if mut_count_id is not None:
        seg_counts = seg_counts[:,mut_count_id]
    seg_mat = gene.segmentgraph.segments[0]
    som_expr_dict = {}

    gidx = countinfo.gene_idx_dict[gene.name]
    seg_pos_list = np.arange(countinfo.gene_id_to_segrange[gidx][0], countinfo.gene_id_to_segrange[gidx][1])
    for ipos in mutation_pos:
        seg_id = bisect.bisect(seg_mat,ipos)
        if seg_id > 0 and ipos <= gene.segmentgraph.segments[1][seg_id-1]: # the mutation is within the pos
            expr = seg_counts[seg_id - 1]
            som_expr_dict[ipos] = expr

    return som_expr_dict


def get_mut_comb(exon_som_dict, vertex_list):
    """
    Get all the mutation combination given the mutation given.
    Parameters
    ----------
    exon_som_dict: dict, keys (exon id), values (somatic mutation position)
    vertex_list: list, array of vertex ids

    Returns
    -------
    mut_comb: list of tuple, list of mutation combination

    """
    mut_comb = [np.nan]
    if exon_som_dict is not None:
        exon_list = map(lambda x: exon_som_dict[x], vertex_list)
        all_comb = get_all_comb(reduce(lambda x, y: x + y, exon_list))
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
    submutation = Mutation(mode=mutation.mode,somatic_mutation_dict=somatic_mutation_sub_dict,germline_mutation_dict=germline_mutation_sub_dict)
    return submutation
