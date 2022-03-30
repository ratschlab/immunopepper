"""Contain functions to deal with mutation"""
import bisect
import logging
import sys
import numpy as np
import pysam

from functools import reduce

from spladder.classes.gene import Gene

from immunopepper.namedtuples import CountInfo, Mutation
from immunopepper.preprocess import parse_mutation_from_maf
from immunopepper.preprocess import parse_mutation_from_vcf
from immunopepper.utils import get_all_comb



def get_mutated_sequence(fasta_file: str, chromosome: str, pos_start: int, pos_end: int,
                         mutation_dict: dict[int:dict[str:str]]):
    """ Reads the specified subsequence from the given fasta file and applies the mutations in mutation_dict.

    :param fasta_file: file containing the reference sequence to be mutated
    :param chrm: the chromosome to which the mutations should be applied
    :param pos_start: starting position for subsequence to be mutated
    :param pos_end: ending position for subsequence to be mutated
    :param mutation_dict: dictionary mapping a position to the mutation variants, e.g.::
        {10 : {'mut_base': '*', 'ref_base': 'A'}}
    :return: the original reference sequence and the mutated sequence in a dictionary. The 'ref' key contains
      the original sequence, and the 'background' key contains the mutated sequence.
    TODO(dd): using a dict here is not right; should be a tuple (or a class)

    """
    with pysam.FastaFile(fasta_file) as fh:
        ref_sequence = fh.fetch(chromosome, pos_start, pos_end)
    return {'ref': ref_sequence,
            'background': _apply_mutations(ref_sequence, pos_start, pos_end, mutation_dict)}


def _apply_mutations(seq: str, pos_start: int, pos_end: int, mut_dict: dict[int:dict[str:str]]):
    """ Apply the given mutations to a reference sequence.  """
    if not mut_dict:
        return seq

    variants = [ipos for ipos in mut_dict.keys() if pos_start <= ipos < pos_end]
    if not variants:
        return seq

    variants.sort()
    # add a sentinel variant at the end of the sequence
    variants.append(pos_start + len(seq))
    mutated_sequences = [seq[:variants[0] - pos_start]]  # copy reference string until the first mutation
    for i in range(len(variants) - 1):  # process all but the last mutation in order
        mut_base = mut_dict[variants[i]]['mut_base']
        ref_base = mut_dict[variants[i]]['ref_base']
        if mut_base != '*':
            mutated_sequences.append(mut_base)
        else:
            mutated_sequences.append(ref_base)
        mutated_sequences.append(seq[variants[i] - pos_start + 1:variants[i + 1] - pos_start])
    return ''.join(mutated_sequences)



def _parse_mutation_file(mutation_tag, mutation_file_path, output_dir, heter_code, target_sample_list=None,
                         mutation_sample=None, name_eq_dict={}):
    """ Reads data from a MAF or VCF file """
    if mutation_file_path.lower().endswith('.maf'):
        return parse_mutation_from_maf(maf_path=mutation_file_path, mutation_mode=mutation_tag,
                                       target_samples=target_sample_list, mutation_sample=mutation_sample,
                                       sample_eq_dict=name_eq_dict, output_dir=output_dir)
    # for VCFs, we also accept hdf5 file format
    elif mutation_file_path.lower().endswith('.vcf') or mutation_file_path.lower().endswith('.h5'):
        return parse_mutation_from_vcf(vcf_path=mutation_file_path, mutation_mode=mutation_tag,
                                       target_samples=target_sample_list, mutation_sample=mutation_sample,
                                       sample_eq_dict=name_eq_dict, heter_code=heter_code,
                                       output_dir=output_dir)

    else:
        logging.error('Unsupported mutation file format: only maf and vcf formats are supported.')
        sys.exit(1)
        
  
def load_mutations(germline_file: str, somatic_file: str, mutation_sample: str, heter_code: int,
                   samples_to_pickle: list[str], sample_name_map_file: str, cache_dir: str):
    """
    Read the given germline/somatic mutation files and return their content as a :class:`Mutation` object.
    :param germline_file: file containing germline mutations, empty if no germline mutations should be applied
    :param somatic_file: file containing somatic mutations, empty if no somatic mutations should be applied
    :param heter_code: Can be 0 or 2. specify which number represents the heterozygous allele (HDF5 VCF only)
        0: 0-> homozygous alternative(1|1), 1-> heterozygous(0|1,1|0) 2->homozygous reference(0|0)
        2: 0-> homozygous reference(0|0), 1-> heterozygous(0|1,1|0) 2->homozygous alternative(1|1)
    :param samples_to_pickle: ids of the samples for which we want to pickle the mutations 
    :param sample_name_map_file: file mapping sample id names in the Spladder graph/count files to
        sample ids in the mutation files. with column 1: graph sample name, 
        column 2: germline sample name, 
        column 3: somatic sample name
    :cache_dir: if not None, pickle result to disk for faster loading
    :return: a :class:`Mutation` instance containing the loaded somatic and germline mutations.
    :rtype Mutation:
    """

    graph_to_somatic_names = {}
    graph_to_germline_names = {}
    if sample_name_map_file is not None:
        with open(sample_name_map_file, "r") as f:
            for line in f.readlines():
                equivalence = line.split("\n")[0].split('\t')
                # Restrict to the samples which should be extracted and saved as pickle
                if samples_to_pickle and (equivalence[0] not in samples_to_pickle):
                    continue
                if len(equivalence) == 2:
                    graph_to_somatic_names[equivalence[0]] = equivalence[1]
                    graph_to_germline_names[equivalence[0]] = equivalence[1]
                elif len(equivalence) == 3:
                    graph_to_germline_names[equivalence[0]] = equivalence[1]
                    graph_to_somatic_names[equivalence[0]] = equivalence[2]
                else:
                    logging.error(f'Invalid sample map file: {sample_name_map_file}.\n\tLine: {line}')
                    sys.exit(1)
    # If no name equivalence file is provided
    else:
        # Take the names of the samples to extract and pickle
        if samples_to_pickle:
            for sample in samples_to_pickle:
                graph_to_somatic_names[sample] = sample
        # Take the names of the mutation sample
        else:
            graph_to_somatic_names[mutation_sample] = mutation_sample
        graph_to_germline_names = graph_to_somatic_names

    mutation_mode = 'ref'
    if germline_file and somatic_file:
        mutation_mode = 'somatic_and_germline'
    elif germline_file:
        mutation_mode = 'germline'
    elif somatic_file:
        mutation_mode = 'somatic'

    somatic_dict = {}
    germline_dict = {}
    if germline_file:
        germline_dict = _parse_mutation_file("germline", germline_file, cache_dir, heter_code,
                                             target_samples, mutation_sample,
                                             graph_to_germline_names)
    if somatic_file:
        somatic_dict = _parse_mutation_file("somatic", somatic_file, cache_dir, heter_code,
                                            target_samples, mutation_sample,
                                            graph_to_somatic_names)

    return Mutation(mutation_mode, germline_dict=germline_dict, somatic_dict=somatic_dict)


def exon_to_mutations(gene: Gene, mutation_pos: dict):
    """
    Builds a dictionary mapping exon ids (vertex indices in the spladder graph) to mutation positions.
    """
    exon_list = gene.splicegraph.vertices
    exon_som_dict = {k: [] for k in range(exon_list.shape[1])}
    for mutation in mutation_pos:
        for i in range(exon_list.shape[1]):
            if mutation in range(exon_list[0, i], exon_list[1, i]):
                exon_som_dict[i].append(mutation)
    exon_som_dict[np.nan] = []  # for single cds case
    return exon_som_dict


def exon_to_expression(gene: Gene, mutation_pos: list[int], count_info: CountInfo, seg_counts: np.ndarray,
                       mut_count_id):
    """
    Builds a dictionary mapping exon ids to expression data.
    """
    if count_info is None:
        return None
    if mut_count_id is not None:
        seg_counts = seg_counts[:, mut_count_id]
    seg_mat = gene.segmentgraph.segments[0]
    som_expr_dict = {}

    for ipos in mutation_pos:
        seg_id = bisect.bisect(seg_mat, ipos)
        if seg_id > 0 and ipos <= gene.segmentgraph.segments[1][seg_id - 1]:  # the mutation is within the pos
            expr = seg_counts[seg_id - 1]
            som_expr_dict[ipos] = expr

    return som_expr_dict


def get_mut_comb(exon_to_mutations: dict[int, list[int]], exon_ids: list[int]):
    """
    Get all the mutation combinations for the given vertices (representing exon ids).
    For example, if exon_to_mutations is {1:[100,200], 2:[300]} and exon_ids is [1,2], the function will return::
         [nan, (100,), (200,), (300,), (100, 200), (100, 300), (200, 300), (100, 200, 300)]
    :param exon_to_mutations: maps exon ids to somatic mutation positions
    :param exon_ids: list of vertices (exon ids) in the spladder graph
    :return: list of tuple, each tuple containing a possible mutation combination.
    """
    mut_comb = [np.nan]
    if exon_to_mutations is not None:
        exon_list = map(lambda x: exon_to_mutations[x], exon_ids)
        all_comb = get_all_comb(reduce(lambda x, y: x + y, exon_list))
        mut_comb += all_comb
    return mut_comb


def get_sub_mutations(mutation: Mutation, sample: str, chromosome: str):
    """ Get the mutations for the given sample and chromosome """
    germline_mutation_sub_dict = mutation.germline_dict.get((sample, chromosome), {})
    somatic_mutation_sub_dict = mutation.somatic_dict.get((sample, chromosome), {})
    return Mutation(mode=mutation.mode, somatic_dict=somatic_mutation_sub_dict,
                    germline_dict=germline_mutation_sub_dict)
