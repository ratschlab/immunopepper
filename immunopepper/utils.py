"""Contain functions to help compute, to preprocess"""
import bisect
import itertools
import logging
import numpy as np
import os
import pandas as pd
import pickle
import psutil
import pyarrow as pa
import signal as sig
import sys


from immunopepper.namedtuples import Idx


def to_adj_list(adj_matrix):
    """
    Converts a binary adjacency matrix to a list of directed edges
    """
    adj_list = []
    assert (adj_matrix.shape[0] == adj_matrix.shape[1])
    for idx in range(adj_matrix.shape[0]):
        for jdx in range(idx, adj_matrix.shape[0]):
            if adj_matrix[idx, jdx] == 1:
                adj_list.append([idx, jdx])
    return adj_list


def get_successor_list(adj_matrix, vertex_map, read_strand):
    """ Returns a list of successors by vertex, sensitive to the read strand

    Parameters
    ----------
    adj_matrix: 2D array with size [len(edges),len(edges)]. '1' represents the two edges
        are connected.
    vertex_map: 2D array with size [2, len(vertexs)]. Vertex positions
    read_strand: The gene strand. '+' or '-'.

    Returns
    -------
    succ_list: List[List[int]]. The connected vertex list for each vertex

    """
    succ_list = []
    assert (adj_matrix.shape[0] == adj_matrix.shape[1])
    for idx in range(adj_matrix.shape[0]):
        succ_list.append([])
        for jdx in range(adj_matrix.shape[0]):
            if adj_matrix[idx, jdx] == 1 and \
               ((read_strand == "+" and vertex_map[0][idx] <= vertex_map[0][jdx]) or \
                (read_strand == "-" and vertex_map[1][idx] >= vertex_map[1][jdx])):
                    succ_list[idx].append(jdx)
    return succ_list


def find_overlapping_cds_simple(v_start, v_stop, cds_begins, strand):
    """
    Find overlapping CDS within an exon given a list of CDS starts
    """
    # cds_start = cds_begin[0]
    if strand == '+':
        return list(filter(lambda x: x[0] >= v_start and x[0] < v_stop, cds_begins))
    else:
        return list(filter(lambda x: x[0] > v_start and x[0] <= v_stop, cds_begins))


def leq_strand(coord1, coord2, strand):
    if strand == "+":
        return coord1 <= coord2
    else:
        return coord1 >= coord2


def encode_chromosome(in_num):
    """  Encodes human chromosome numbers to strings (23==X, 24==Y, 25 ==MT, the rest remain unchanged. """
    convert_dict = {23: "X", 24: "Y", 25: "MT"}
    return convert_dict[in_num] if in_num in convert_dict else str(in_num)


def get_sub_mut_dna(background_seq, coord, variant_comb, somatic_mutation_sub_dict, strand, gene_start):
    """ Get the mutated dna sub-sequence according to mutation specified by the variant_comb.

    Parameters
    ----------
    background_seq: List(str). backgound sequence.
    coord: Namedtuple containing the vertex coordinates.
    variant_comb: List(int). List of variant position. Like ['38', '43']
    somatic_mutation_sub_dict: Dict. variant position -> variant details.
    strand: gene strand

    Returnvariant_combs
    -------
    sub_dna: str. dna when applied somatic mutation.

    """
    def _get_variant_pos_offset(variant_pos, coord_pair_list, strand):
        offset = 0
        takes_effect = False
        for p1,p2 in coord_pair_list:
            if variant_pos >= p1 and variant_pos < p2:
                if strand == '+':
                    offset += variant_pos - p1
                else:
                    offset += p2 - variant_pos - 1
                takes_effect = True
                break
            else:
                offset = p2 - p1

        return offset if takes_effect else np.nan

    real_coord = list(filter(lambda x: x is not np.nan and x != None, coord))
    assert len(real_coord) % 2 == 0
    coord_pair_list = list(zip(real_coord[::2], real_coord[1::2]))

    if strand == '+':
        sub_dna = ''.join([background_seq[pair[0] - gene_start:pair[1] - gene_start] for pair in coord_pair_list])
    else:
        sub_dna = ''.join([background_seq[pair[0] - gene_start:pair[1] - gene_start][::-1] for pair in coord_pair_list])

    if variant_comb is np.nan : # no mutation exist
        return sub_dna

    relative_variant_pos = [_get_variant_pos_offset(variant_ipos, coord_pair_list, strand) for variant_ipos in variant_comb]
    for i,variant_ipos in enumerate(variant_comb):
        mut_base = somatic_mutation_sub_dict[variant_ipos]['mut_base']
        ref_base = somatic_mutation_sub_dict[variant_ipos]['ref_base']
        pos = relative_variant_pos[i]
        if pos is not np.nan:
            sub_dna = sub_dna[:pos] + mut_base + sub_dna[pos+1:]
    return sub_dna


def get_size_factor(samples, lib_file_path):
    """ Get the expression adjustment parameter for certain samples.

    Parameters
    ----------
    samples: List[str]. List of samples.
    lib_file_path: str. libsize file path

    Returns
    -------
    sf: 1darray. size factor.

    """
    libs = np.loadtxt(lib_file_path, dtype='str', skiprows=1, delimiter='\t')
    a, b = np.where(samples[:, np.newaxis] == libs[:, 0])
    assert np.all(libs[b, 0] == samples)
    libs = libs[b, :]
    med = np.median(libs[:, 1].astype('float'))
    sf = med / libs[:, 1].astype('float')
    return sf


def get_all_comb(array, r=None):
    """ Return all subsequences of the given array with length smaller than r.
    For example, get_all_comb([1,2,3]) will return::
       [(1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)],
    and get_all_comb([1,2,3], 2) will return:
        [(1,), (2,), (3,), (1, 2), (1, 3), (2, 3)]
    :param array: array to get subsequences from
    :param r: the maximum length of returned subsequences, or None for no limit
    :return: list of subsequences
    :rtype: list[tuple]
    """
    if r is None:
        r = len(array)
    return [_ for i in range(1, r + 1) for _ in itertools.combinations(array, i)]


def is_mutation_within_exon_pair(pos_start, pos_end, mutation_pos):
    return list(filter(lambda x: x > pos_start and x < pos_end, mutation_pos))


def is_isolated_cds(gene, gene_info, idx):
    """
    Indicate if a peptide is translated by an isolated cds

    Currently, it is a simple version. If the exon where the cds is
    has no connected exon, it is isolated. (the method is not true in some cases)

    Parameters
    ----------
    gene: gene object
    idx: exon id

    """

    if len(gene_info.vertex_succ_list[idx]) > 0:
        return False

    return np.sum(gene.splicegraph.edges[:, idx]) == 0


def get_exon_expr(gene, vstart, vstop, countinfo, Idx, seg_counts):
    """ Split the exon into segments and get the corresponding counts.

    Parameters
    ----------
    gene: Object. Generated by SplAdder.
    vstart: int. Start position of given exon.
    vstop: int. Stop position of given exon.
    countinfo: Namedtuple, containing spladder count information
    Idx: Namedtuple, has attribute idx.gene and idx.sample
    seg_counts: np.array, array of spladder segment counts from gene and sample of interest

    Returns
    -------
    expr_list: List[Tuple(int,float)] (int, float) represents the length of segment
        and the expression count of that segment.

    """
    if countinfo is None:
        return np.zeros((0, 1), dtype='float') #[np.nan]

    out_shape = (seg_counts.shape[1] + 1) if len(seg_counts.shape) > 1 else 2
    if vstart is np.nan or vstop is np.nan:  # isolated exon case
        return np.zeros((0, out_shape), dtype='float')


    segments = gene.segmentgraph.segments

    sv1_id = bisect.bisect(segments[0], vstart) - 1
    sv2_id = bisect.bisect(segments[0], vstop) - 1
    if sv1_id == sv2_id:
        if len(seg_counts.shape) > 1:
            expr_list = np.c_[np.array([vstop - vstart]), [seg_counts[sv1_id, :]]]
        else:
            expr_list = np.array([(vstop - vstart, seg_counts[sv1_id])])
    else:
        if len(seg_counts.shape) > 1:
            expr_list = np.c_[segments[1, sv1_id:sv2_id + 1] - segments[0, sv1_id:sv2_id + 1], seg_counts[sv1_id:sv2_id + 1, :]]
        else:
            expr_list = np.c_[segments[1, sv1_id:sv2_id + 1] - segments[0, sv1_id:sv2_id + 1], seg_counts[sv1_id:sv2_id + 1, np.newaxis]]
        expr_list[0, 0] -= (vstart - segments[0, sv1_id])
        expr_list[-1, 0] -= (segments[1, sv2_id] - vstop)
        if gene.strand == '-': # need to reverse epression list to match the order of translation
            expr_list = expr_list[::-1]
    return expr_list

def get_segment_expr(gene, coord, countinfo, Idx, seg_counts, cross_graph_expr):
    """ Get the segment expression for one exon-pair.
    Apply 'get_exon_expr' for each exon and concatenate them.

    Parameters
    ----------
    gene: Object. Generated by SplAdder.
    coord: NamedTuple. has attribute ['start_v1', 'stop_v1', 'start_v2', 'stop_v2']
        contains the true four position of exon pairs (after considering read framee)
        that outputs the peptide.
    countinfo: Namedtuple, contains count information generated by SplAdder
        has attributes [sample_idx_dict, gene_idx_dict, gene_id_to_segrange, gene_id_to_edgerange, h5fname]
    Idx: Namedtuple, has attribute idx.gene and idx.sample
    seg_counts: np.array, array of spladder segment counts from gene and sample of interest


    Returns
    -------
    expr_meta_file: float. The average expression counts for the given exon-pair.
        Is nan if matrix mode, because we do not report peptide expressions per sample in this mode
    expr1: List[Tuple(int,float)] (int, float) represents the length of segment
        and the expression count of that segment.

    """
    expr_meta_file = np.nan
    if coord.start_v3 is None:
        expr_list = np.vstack([get_exon_expr(gene, coord.start_v1, coord.stop_v1, countinfo, Idx, seg_counts ),
                               get_exon_expr(gene, coord.start_v2, coord.stop_v2, countinfo, Idx, seg_counts )])
    else:
        expr_list = np.vstack([get_exon_expr(gene, coord.start_v1, coord.stop_v1, countinfo, Idx, seg_counts ),
                               get_exon_expr(gene, coord.start_v2, coord.stop_v2, countinfo, Idx, seg_counts ),
                               get_exon_expr(gene, coord.start_v3, coord.stop_v3, countinfo, Idx, seg_counts )])

    seg_len = np.sum(expr_list[:, 0])

    n_samples = expr_list[:, 1:].shape[1]
    len_factor = np.tile(expr_list[:, 0], n_samples).reshape(n_samples, expr_list.shape[0]).transpose()
    mean_expr = (np.sum(expr_list[:, 1:]*len_factor, 0) / seg_len).astype(int) if seg_len > 0 else np.zeros(n_samples).astype(int)
    if not cross_graph_expr:
        expr_meta_file = mean_expr[0]
    return expr_meta_file, expr_list


def get_total_gene_expr(gene, countinfo, Idx, seg_expr, cross_graph_expr):
    """ get total reads count for the given sample and the given gene
    actually total_expr = reads_length*total_reads_counts
    """
    if len(seg_expr.shape) == 1:
        n_samples = 1
    else:
        n_samples = seg_expr.shape[1]

    if countinfo is None:
        return [np.nan] * n_samples
    seg_len = gene.segmentgraph.segments[1] - gene.segmentgraph.segments[0]

    if cross_graph_expr:
        total_expr = np.sum(seg_len * seg_expr.T, axis=1)
        total_expr = total_expr.tolist()
    else:
        total_expr = [np.sum(seg_len*seg_expr)]
    return total_expr


def get_idx(countinfo, sample, gene_idx):
    """ Create a  aggregated Index with namedtuple idx
    Combine the gene_idx, sample_idx

    Parameters
    ----------
    countinfo: NameTuple, containing spladder count information
    sample: str. sample name
    gene_idx: int. Gene index, mainly for formatting the output.

    """
    if not countinfo is None:
        if not sample in countinfo.sample_idx_dict:
            sample_idx = None
            if sample != 'cohort':
                logging.warning("utils.py: The sample {} is not in the count file. Program proceeds without outputting expression data.".format(sample))
        else:
            sample_idx = countinfo.sample_idx_dict[sample]
    else:
        sample_idx = None

    return Idx(gene_idx, sample_idx)


def create_libsize(expr_distr_fp, output_fp, sample, debug=False):
    """ create library_size text file.

    Calculate the 75% expression and sum of expression for each sample
    and write into output_fp.

    Parameters
    ----------
    expr_distr_dict: Dict. str -> List(float). Mapping sample to the expression of all exon pairs
    output_fp: file pointer. library_size text
    debug: Bool. In debug mode, return the libsize_count dictionary.
    """
    libsize_exists = ''
    sample_expr_distr = sample_expr_distr = pd.read_csv(expr_distr_fp['path'], sep='\t', compression='gzip') #TODO gather all the partitions here

    # change types
    convert_dict = {}
    for col in sample_expr_distr.columns:
        if col != 'gene':
            convert_dict[col] = float
    sample_expr_distr = sample_expr_distr.astype(convert_dict)

    # compute libsizes
    df_libsize = pd.DataFrame({'sample': sample_expr_distr.columns[1:],
                                 'libsize_75percent': np.percentile(sample_expr_distr.iloc[:, 1:], 75, axis=0, interpolation='linear'),
                                  'libsize_total_count': np.sum(sample_expr_distr.iloc[:, 1:], axis=0)}, index = None)


    if os.path.isfile(output_fp):
        previous_libsize = pd.read_csv(output_fp, sep = '\t')
        df_libsize = pd.concat([previous_libsize, df_libsize], axis=0).drop_duplicates(subset=['sample'], keep='last')
        libsize_exists = ': append to existing file.'
    logging.info(f'Saved library size results to {output_fp}{libsize_exists}')
    df_libsize.to_csv(output_fp, sep='\t', index=False)


def get_concat_peptide(front_coord_pair, back_coord_pair, front_peptide, back_peptide, strand, k=None):
    """
    Get the concatenated peptide from possible match peptide.
    Parameters
    ----------
    front_coord_pair: str. coordinate string for the front exon pair.
    back_coord_pair: str. coordinate string for the back exon pair.
    front_peptide: str. peptide translated from the front exon pair.
    back_peptide: str. peptide translated from the back exon pair.
    strand: str. '+' or '-'
    k: k for k-mer.

    Returns
    -------
    new_peptide: str. concatenated peptide. If none is found, return empty string

    Examples
    --------
    front_pep: 'MKTT', back_pep: 'TTAC', concatenated_pep: 'MKTTAC'

    """
    def get_longest_match_position(front_str,back_str,L=None):
        if L is None:
            L = min(len(front_str),len(back_str))
        for i in reversed(list(range(1,L+1))):
            if front_str[-i:] == back_str[:i]:
                return i
        return None
    if strand == '+':
        front_coord = front_coord_pair.stop_v2
        back_coord = back_coord_pair.start_v1
    else:
        front_coord = front_coord_pair.start_v2
        back_coord = back_coord_pair.stop_v1
    if abs(front_coord-back_coord) % 3 == 0:
        if front_coord == back_coord:  # no intersection and we concatenate them directly
            new_peptide = front_peptide + back_peptide
        else:
            pep_common_num = get_longest_match_position(front_peptide,back_peptide,L=k)
            if pep_common_num is None:
                new_peptide = ''
            else:
                new_peptide = front_peptide + back_peptide[pep_common_num:]
        return new_peptide
    else:
        return ''

def replace_I_with_L(kmer):
    return kmer.replace('I', 'L')

def check_chr_consistence(ann_chr_set, mutation, graph_data):
    germline_chr_set = set()
    somatic_chr_set = set()
    mode = mutation.mode
    if mutation.germline_dict:
        germline_chr_set = set([item[1] for item in mutation.germline_dict.keys()])
    if mutation.somatic_dict:
        somatic_chr_set = set([item[1] for item in mutation.somatic_dict.keys()])
    whole_mut_set = somatic_chr_set.union(germline_chr_set)
    common_chr = whole_mut_set.intersection(ann_chr_set)

    if mode != 'ref' and len(common_chr) == 0:
        logging.error("Mutation file has different chromosome naming from annotation file, please check")
        sys.exit(0)
    if len(graph_data) > 0:
        gene_chr_set = set([gene.chr for gene in graph_data])
        new_chr_set = gene_chr_set.difference(ann_chr_set)
        if len(new_chr_set) > 0:
            logging.error("Gene object has different chromosome naming from annotation file, please check")
            sys.exit(0)


def unpickler(picklefile):
    try:
        while True:
            yield pickle.load(picklefile)
    except EOFError:
        pass


def print_memory_diags(disable_print=False):
    """
    Print memory diagnostics including the active resident set size
    """
    process = psutil.Process(os.getpid())
    memory = process.memory_info().rss/1000000000.0
    if not disable_print:
        logging.info('\tMemory usage: {:.3f} GB'.format(memory))
    return memory

def pool_initializer():
    return sig.signal(sig.SIGINT, sig.SIG_IGN)
