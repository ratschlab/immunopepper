import gzip
import os
import pickle

import Bio.SeqIO as BioIO
import pytest

from immunopepper.immuno_mutation import apply_germline_mutation,mut_replace
from immunopepper.immuno_preprocess import preprocess_ann, genes_preprocess, \
    parse_mutation_from_vcf, parse_mutation_from_maf
from immunopepper.utils import get_sub_mut_dna,complementary_seq
from immunopepper.io_utils import load_pickled_graph
from immunopepper.immuno_filter import get_true_variant_comb
from immunopepper.immuno_model import create_output_kmer

from collections import namedtuple
data_dir = os.path.join(os.path.dirname(__file__), 'test_base','data')

@pytest.fixture
def load_gene_data():
    f = open(os.path.join(data_dir, 'posgraph','spladder',
                          'genes_graph_conf3.merge_graphs.pickle'), 'r')
    ann_path = os.path.join(data_dir, 'test_basepos.gtf')
    ref_path = os.path.join(data_dir, 'test_basepos.fa')

    (graph_data, graph_meta) = load_pickled_graph(f)  # cPickle.load(f)
    genetable = preprocess_ann(ann_path)
    interesting_chr = map(str, range(1, 23)) + ["X", "Y", "MT"]
    seq_dict = {}
    for record in BioIO.parse(ref_path, "fasta"):
        if record.id in interesting_chr:
            seq_dict[record.id] = str(record.seq).strip()

    gene = graph_data[0]
    chrm = gene.chr.strip()
    ref_seq = seq_dict[chrm]
    return graph_data, ref_seq, genetable.gene_to_cds_begin


@pytest.fixture
def load_mutation_data(testname):
    data_dir = os.path.join(os.path.dirname(__file__), testname, 'data')
    vcf_filename = testname + 'pos.vcf'
    maf_filename = testname + 'pos.maf'
    vcf_path = os.path.join(data_dir, vcf_filename)
    maf_path = os.path.join(data_dir, maf_filename)
    mutation_dic_vcf = parse_mutation_from_vcf(vcf_path)
    mutation_dic_maf = parse_mutation_from_maf(maf_path)
    return mutation_dic_vcf, mutation_dic_maf


@pytest.fixture
def load_ref_seq():
    pos_ref_seq = 'GTAATGTGTAAGATGACGCACGCATGGTGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTTCAG' \
              'GTACGTATAGTTGCTAAGTAGTCGACGTTCTGGTGGTAGACCTTTTGTTACCCAATTGGGTAAGGCAGCCATGTTGATATACAT'
    neg_ref_seq = 'ATGTATATCAACATGGCTGCCTTACCCAATTGGGTAACAAAAGGTCTACCACCAGAACGTCGACTAC' \
                  'TTAGCAACTATACGTACCTGAACTCGAACTTACTCCGCAACCCATCTCCAATACCACCATGCGTGCGTCATCTTACACATTAC'
    return pos_ref_seq,neg_ref_seq

@pytest.fixture
def load_mut_dict():

    pos_var_del = {'ref_base': 'CA', 'mut_base': 'C', 'strand': '+',
                'Variant_Classification': 'Silent', 'Variant_Type': 'DEL'}
    pos_var_snp = {'ref_base': 'A', 'mut_base': 'G', 'strand': '+',
                'Variant_Classification': 'Silent', 'Variant_Type': 'SNP'}
    pos_var_ins = {'ref_base': 'G', 'mut_base': 'GT', 'strand': '+',
                'Variant_Classification': 'Silent', 'Variant_Type': 'INS'}
    pos_mutation_dic = {}
    pos_mutation_dic[18] = pos_var_del
    pos_mutation_dic[38] = pos_var_snp
    pos_mutation_dic[45] = pos_var_ins

    neg_var_del = {'ref_base': 'TG', 'mut_base': 'G', 'strand': '-',
                'Variant_Classification': 'Silent', 'Variant_Type': 'DEL'}
    neg_var_snp = {'ref_base': 'T', 'mut_base': 'C', 'strand': '-',
                'Variant_Classification': 'Silent', 'Variant_Type': 'SNP'}
    neg_var_ins = {'ref_base': 'C', 'mut_base': 'AC', 'strand': '-',
                'Variant_Classification': 'Silent', 'Variant_Type': 'INS'}

    neg_mutation_dic = {}
    neg_mutation_dic[104] = neg_var_ins
    neg_mutation_dic[111] = neg_var_snp
    neg_mutation_dic[130] = neg_var_del
    return pos_mutation_dic, neg_mutation_dic


def test_preprocess(load_gene_data):
    graph_data, seq_dict, gene_cds_begin_dict = load_gene_data
    genes_preprocess(graph_data, gene_cds_begin_dict)
    assert graph_data[0].nvertices == 8


def test_get_sub_mut_dna(load_ref_seq,load_mut_dict):
    # --------------
    # positive case
    pos_ref_seq, neg_ref_seq = load_ref_seq
    pos_mut_dict, neg_mut_dict = load_mut_dict


    pos_groundtruth = ['GATGACGCCGCATGGTGATGGGTTGCGGA',
                       'GATGACGCCGCATGGTGGTGGGTTGCGGA',
                       'GATGACGCCGCATGGTGATGGGTTGTCGGA',
                       'GATGACGCCGCATGGTGGTGGGTTGTCGGA']
    vlist = [11, 29, 38, 50]
    variant_comb_list = [(18,), (18, 38), (18, 45), (18, 38, 45)]
    for i, variant_comb in enumerate(variant_comb_list):
        sub_dna = get_sub_mut_dna(pos_ref_seq, vlist[0], vlist[1], vlist[2],
                                  vlist[3], variant_comb,
                                  pos_mut_dict, '+')
        assert sub_dna == pos_groundtruth[i]

    # --------------
    # negative case
    vlist = [121, 139, 100, 112]
    variant_comb_list = [(130,), (130, 111), (111, 104), (130, 111, 104)]

    # neg_dna_list = list(neg_ref_seq)
    # assert neg_dna_list[111] == 'T'
    # neg_dna_list[111] = 'C'
    # assert ''.join(neg_dna_list[104]) == 'C'
    # neg_dna_list[104] = 'AC'
    # assert ''.join(neg_dna_list[130:132]) == 'TG'
    # neg_dna_list[130:132] = ['','G']
    # complementary_seq(''.join(neg_dna_list[121:139][::-1]+neg_dna_list[100:112][::-1]))

    neg_groundtruth = ['GATGACGCCGCATGGTGATGGGTTGCGGA',
                       'GATGACGCCGCATGGTGGTGGGTTGCGGA',
                       'GATGACGCACGCATGGTGGTGGGTTTGCGGA',
                       'GATGACGCCGCATGGTGGTGGGTTTGCGGA']
    for i, variant_comb in enumerate(variant_comb_list):
        sub_dna = complementary_seq(get_sub_mut_dna(neg_ref_seq, vlist[0], vlist[1], vlist[2],
                                  vlist[3], variant_comb,
                                  neg_mut_dict, '-'))
        assert sub_dna == neg_groundtruth[i]


def test_mut_replace():
    ori_seq = 'CGTACAATGA'
    variant_pos = 2
    ref_base = 'TA'
    mut_base = 'T'
    mut_seq = mut_replace(ori_seq,variant_pos,ref_base,mut_base)
    assert mut_seq == ['C', 'G', 'T', '', 'C', 'A', 'A', 'T', 'G', 'A']

    variant_pos = 2
    ref_base = 'T'
    mut_base = 'TA'
    mut_seq = mut_replace(ori_seq,variant_pos,ref_base,mut_base)
    assert mut_seq == ['C', 'G', 'TA', 'A', 'C', 'A', 'A', 'T', 'G', 'A']


def test_apply_germline_mutation(load_ref_seq,load_mut_dict):
    pos_ref_seq, neg_ref_seq = load_ref_seq
    pos_mut_dict, neg_mut_dict = load_mut_dict
    ref_mut_seq = apply_germline_mutation(pos_ref_seq,0,len(pos_ref_seq),pos_mut_dict)
    pos_dna = 'GTAATGTGTAAGATGACGCACGCATGGTGGTATTGGAGATGGGT' \
              'TGCGGAGTAAGTTCGAGTTCAGGTACGTATAGTTGCTAAGTAGTCGACGTTCTGGTGGTAGACCTTTTGTTACCCAATTGGGTAAGGCAGCCATGTTGATATACAT'
    pos_dna_list = list(pos_dna)
    assert pos_dna_list[18:20] == ['C','A']
    pos_dna_list[18:20] = ['C','']
    assert pos_dna_list[38] == 'A'
    pos_dna_list[38] = 'G'
    assert pos_dna_list[45] == 'G'
    pos_dna_list[45] = 'GT'
    assert ref_mut_seq['ref'] == list(pos_ref_seq)
    assert ref_mut_seq['background'] == pos_dna_list

    neg_dna = 'ATGTATATCAACATGGCTGCCTTACCCAATTGGGTAACAAAAGGTCTACCACCAGAACGTCGACTACTTAGCAACTATACGTACCTGAACTCGAACTTACTCCG' \
              'CAACCCATCTCCAATACCACCATGCGTGCGTCATCTTACACATTAC'
    neg_dna_list = list(neg_dna)
    assert neg_dna_list[130:132] == ['T','G']
    neg_dna_list[130:132] = ['G','']
    assert neg_dna_list[111] == 'T'
    neg_dna_list[111] = 'C'
    assert neg_dna_list[104] == 'C'
    neg_dna_list[104] = 'AC'
    ref_mut_seq = apply_germline_mutation(neg_ref_seq,0,len(neg_ref_seq),neg_mut_dict)
    assert ref_mut_seq['ref'] == list(neg_ref_seq)
    assert ref_mut_seq['background'] == neg_dna_list


def test_complementary_seq():
    dna = 'AGCT'
    assert complementary_seq(dna) == 'TCGA'
    assert complementary_seq(list(dna)) == 'TCGA'
    false_dna = 'aGcT'
    assert complementary_seq(false_dna) == 'aCcA'


def test_get_true_variant_comb():
    variant_comb = '38;43'.split(';')
    exon_coord = '39;50;66;74'.split(';')
    true_variant = get_true_variant_comb(variant_comb,exon_coord)
    assert true_variant == '43'
    variant_comb = '106;111'.split(';')
    exon_coord = '100;110;76;84'.split(';')
    true_variant = get_true_variant_comb(variant_comb,exon_coord)
    assert true_variant == '106'
    variant_comb = '.'.split(';')
    exon_coord = '100;110;76;84'.split(';')
    true_variant = get_true_variant_comb(variant_comb,exon_coord)
    assert true_variant == '.'
