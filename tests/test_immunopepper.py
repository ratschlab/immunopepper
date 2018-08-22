import gzip
import os
import pickle

import Bio.SeqIO as BioIO
import pytest

from immunopepper.immuno_mutation import apply_germline_mutation
from immunopepper.immuno_preprocess import preprocess_ann, genes_preprocess, \
    parse_mutation_from_vcf, parse_mutation_from_maf
from immunopepper.utils import get_sub_mut_dna
from immunopepper.io_utils import load_pickled_graph

data_dir = os.path.join(os.path.dirname(__file__), 'data')


@pytest.fixture
def load_gene_data():
    f = open(os.path.join(data_dir, 'spladder',
                          'genes_graph_conf3.merge_graphs.pickle'), 'r')
    ann_path = os.path.join(data_dir, 'test1.gtf')
    ref_path = os.path.join(data_dir, 'test1.fa')

    (graph_data, graph_meta) = load_pickled_graph(f)  # cPickle.load(f)
    gene_cds_begin_dict, gene_to_transcript_table, transcript_to_cds_table = preprocess_ann(
        ann_path)
    interesting_chr = map(str, range(1, 23)) + ["X", "Y", "MT"]
    seq_dict = {}
    for record in BioIO.parse(ref_path, "fasta"):
        if record.id in interesting_chr:
            seq_dict[record.id] = str(record.seq).strip()

    gene = graph_data[0]
    chrm = gene.chr.strip()
    ref_seq = seq_dict[chrm]
    return graph_data, ref_seq, gene_cds_begin_dict


@pytest.fixture
def load_mutation_data():
    vcf_path = os.path.join(data_dir, 'test1.vcf')
    maf_path = os.path.join(data_dir, 'test1.maf')
    mutation_dic_vcf = parse_mutation_from_vcf(vcf_path)
    mutation_dic_maf = parse_mutation_from_maf(maf_path)

    return mutation_dic_vcf, mutation_dic_maf


def test_preprocess(load_gene_data):
    graph_data, seq_dict, gene_cds_begin_dict = load_gene_data
    genes_preprocess(graph_data, gene_cds_begin_dict)
    assert graph_data[0].nvertices == 6


def test_germline_mutation(load_gene_data, load_mutation_data):
    graph_data, ref_seq, gene_cds_begin_dict = load_gene_data
    mutation_dic_vcf, mutation_dic_maf = load_mutation_data
    gene = graph_data[0]
    mutation_sub_dic_vcf = mutation_dic_vcf['test1', gene.chr]
    ref_mut_seq = apply_germline_mutation(ref_sequence=ref_seq,
                                          pos_start=gene.start,
                                          pos_end=gene.stop,
                                          mutation_sub_dic_vcf=mutation_sub_dic_vcf)
    assert 'ref' in ref_mut_seq.keys()


def test_get_sub_mut_dna(load_gene_data, load_mutation_data):
    graph_data, ref_seq, gene_cds_begin_dict = load_gene_data
    mutation_dic_vcf, mutation_dic_maf = load_mutation_data
    gene = graph_data[0]

    # made modification to the mutation_dic_vcf
    var_dict = {'ref_base': 'G', 'mut_base': 'A', 'strand': '+',
                'Variant_Classification': 'Silent', 'Variant_Type': 'SNP'}
    mutation_sub_dic_maf = mutation_dic_maf['test1', gene.chr]
    mutation_sub_dic_maf[41] = var_dict
    test_list = [[11, 29, 38, 50],
                 [11, 29, 38, 50],
                 [60, 75, 87, 102],
                 [274, 286, 250, 262],
                 [225, 240, 198, 213]]
    groundtruth = ['GATGACGCACGCATACAGATAGGTAGCGGA',
                   'GATGACGCACGCATACAGATAAGTAGCGGA',
                   'CTTCGGGTACGTATATCGACGTTAGGGTGG',
                   'CTGCGTGCGTATTATCCATCGCCT',
                   'GAAGCCCATGCATATAGCTGCAATCCCACC'
                   ]
    variant_comb = [(40,), (40, 41), '.', (259,), '.']
    strand = ['+', '+', '+', '-', '-']
    for i, vlist in enumerate(test_list):
        sub_dna = get_sub_mut_dna(ref_seq, vlist[0], vlist[1], vlist[2],
                                  vlist[3], variant_comb[i],
                                  mutation_sub_dic_maf, strand[i])
        assert sub_dna == groundtruth[i]
