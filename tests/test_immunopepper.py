import gzip
import os
import pickle

import Bio.SeqIO as BioIO
import pytest

from immunopepper.immuno_mutation import apply_germline_mutation,construct_mut_seq_with_str_concat,get_mutation_mode_from_parser
from immunopepper.immuno_preprocess import preprocess_ann, genes_preprocess, \
    parse_mutation_from_vcf, parse_mutation_from_maf
from immunopepper.utils import get_sub_mut_dna,get_concat_peptide
from immunopepper.io_utils import load_pickled_graph
from immunopepper.main_immuno import parse_arguments
from immunopepper.immuno_model import create_output_kmer

data_dir = os.path.join(os.path.dirname(__file__), 'test1','data')


@pytest.fixture
def load_gene_data():
    f = open(os.path.join(data_dir, 'posgraph','spladder',
                          'genes_graph_conf3.merge_graphs.pickle'), 'r')
    ann_path = os.path.join(data_dir, 'test1pos.gtf')
    ref_path = os.path.join(data_dir, 'test1pos.fa')

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
def load_mutation_data():
    vcf_path = os.path.join(data_dir, 'test1pos.vcf')
    maf_path = os.path.join(data_dir, 'test1pos.maf')
    mutation_dic_vcf = parse_mutation_from_vcf(vcf_path,['test1pos'])
    mutation_dic_maf = parse_mutation_from_maf(maf_path)

    return mutation_dic_vcf, mutation_dic_maf


def test_preprocess(load_gene_data):
    graph_data, seq_dict, gene_cds_begin_dict = load_gene_data
    genes_preprocess(graph_data, gene_cds_begin_dict)
    assert graph_data[0].nvertices == 8


def test_germline_mutation(load_gene_data, load_mutation_data):
    graph_data, ref_seq, gene_cds_begin_dict = load_gene_data
    mutation_dic_vcf, mutation_dic_maf = load_mutation_data
    gene = graph_data[0]
    mutation_sub_dic_vcf = mutation_dic_vcf['test1pos', gene.chr]
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
    mutation_sub_dic_maf = mutation_dic_maf['test1pos', gene.chr]
    mutation_sub_dic_maf[41] = var_dict
    test_list = [[11, 29, 38, 50],
                 [11, 29, 38, 50],
                 [60, 75, 87, 102]]
    groundtruth = ['GATGACGCACGCATGGTGGTGGGTTGCGGA',
                   'GATGACGCACGCATGGTGGTGAGTTGCGGA',
                   'GTTCAGGTACGTATATCGACGTTCTGGTGG',
                   ]
    variant_comb = [(38,), (38, 41), '.']
    strand = ['+', '+', '+']
    for i, vlist in enumerate(test_list):
        sub_dna = get_sub_mut_dna(ref_seq, vlist[0], vlist[1], vlist[2],
                                  vlist[3], variant_comb[i],
                                  mutation_sub_dic_maf, strand[i])
        assert sub_dna == groundtruth[i]


def test_reading_gtf_and_gff3_file():
    gff3_path = os.path.join(data_dir,'small.gencode.v29.gff3')
    gtf_path = os.path.join(data_dir, 'small.gencode.v29.gtf')
    gene_table_gtf = preprocess_ann(gtf_path)
    gene_table_gff = preprocess_ann(gff3_path)
    assert gene_table_gff.gene_to_ts == gene_table_gtf.gene_to_ts
    assert gene_table_gtf.ts_to_cds.keys() == gene_table_gff.ts_to_cds.keys()
    assert gene_table_gtf.ts_to_cds['ENST00000335137.4'] == [(69090, 70005, 0)]
    assert gene_table_gff.ts_to_cds['ENST00000335137.4'] == [(69090, 70008, 0)] # include stop codon in gff3


def test_reading_gtf_and_gff_file():
    gff_path = os.path.join(data_dir,'small.gencode.v19.gff')
    gtf_path = os.path.join(data_dir, 'small.gencode.v19.gtf')
    gene_table_gtf = preprocess_ann(gtf_path)
    gene_table_gff = preprocess_ann(gff_path)
    assert gene_table_gff.gene_to_ts == gene_table_gtf.gene_to_ts
    assert gene_table_gtf.ts_to_cds.keys() == gene_table_gff.ts_to_cds.keys()
    assert gene_table_gff.ts_to_cds['ENST00000335137.3'] == [(69090, 70006, 0)]
    assert gene_table_gtf.ts_to_cds['ENST00000335137.3'] == [(69090, 70005, 0)]


def test_reading_vcf_h5():
    vcf_dict_default_heter_code0 = parse_mutation_from_vcf(os.path.join(data_dir,'test1vcf.h5'),['test1pos','test1neg'])
    vcf_dict_heter_code2 = parse_mutation_from_vcf(os.path.join(data_dir,'test1vcf.h5'),['test1pos','test1neg'],heter_code=2)
    assert len(vcf_dict_default_heter_code0) == 2
    assert vcf_dict_default_heter_code0['test1neg', 'X'][135] == {'mut_base': 'G', 'ref_base': 'C'}
    assert vcf_dict_default_heter_code0['test1pos', 'X'][14] == {'mut_base': 'C', 'ref_base':'G'}
    assert vcf_dict_heter_code2['test1pos', 'X'][135] == {'mut_base': 'G', 'ref_base': 'C'}
    assert vcf_dict_heter_code2['test1neg', 'X'][14] == {'mut_base': 'C', 'ref_base':'G'}
    assert vcf_dict_heter_code2['test1neg', 'X'][135] == {'mut_base': 'G', 'ref_base':'C'}


def test_construct_mut_seq_with_str_concat():
    ref_seq = 'GTAATGTGTAAGATGACGCACGCATGGTGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTTC'
    gt_mut_seq2 = 'GTAATGTGTAAGATGACGCACGCATA'+'C'+'TGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTTC'
    gt_mut_seq3 = 'GTAATGTGTAAGATGACGCACGCATA'+'G'+'TGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTTC'
    mut_dict = {}
    mut_dict[10] = {'mut_base':'*','ref_base':'A'}
    mut_dict[25] = {'mut_base':'A','ref_base':'G'}
    mut_dict[26] = {'mut_base':'C','ref_base':'G'}
    mut_seq1 = construct_mut_seq_with_str_concat(ref_seq, 45, 46, mut_dict) # test unclear mut_base
    assert mut_seq1 == ref_seq
    mut_seq2 = construct_mut_seq_with_str_concat(ref_seq, 25, 27, mut_dict) # [25,27) include 26
    assert mut_seq2 == gt_mut_seq2
    mut_seq3 = construct_mut_seq_with_str_concat(ref_seq, 25, 26, mut_dict) # [25,26) not include 26
    assert mut_seq3 == gt_mut_seq3



def test_get_mutation_mode_from_parser():
    my_args1 = ['--vcf_path', os.path.join(data_dir,'test1pos.vcf'),
               '--maf_path', os.path.join(data_dir,'test1pos.maf'),
               '--mutation_mode', 'somantic']  # bad mutation mode
    args = parse_arguments(my_args1)
    try:
        get_mutation_mode_from_parser(args)
    except SystemExit:
        assert 1
    my_args2 = ['--vcf_path', os.path.join(data_dir,'test1pos.vcf'),
               '--mutation_mode', 'somatic']  # mismatch mutation mode and input files
    args = parse_arguments(my_args2)
    try:
        get_mutation_mode_from_parser(args)
    except SystemExit:
        assert 1


def test_create_output_kmer():
    peptide_list = ['1\nRTHDGLRSTYI','2\nMTHAW']
    expr_lists = [[(8,1000),(1,220),(0,0)]]
    k = 3
    try:
        c = create_output_kmer(peptide_list,expr_lists,k)
    except AssertionError: # make sure len(peptide) == len(expr_lists)
        assert 1
    peptide_list = ['1\nMTHAW']
    expr_lists = [[(8,1000),(1,220),(6,0)]] # test 0 expression
    c = create_output_kmer(peptide_list, expr_lists, k)
    true_output = ['MTH\t1\t913.33', 'THA\t1\t580.0', 'HAW\t1\t246.67']
    assert c == true_output
    expr_lists = [[(8,1000),(1,220),(0,0)]] # test 0 expression
    c = create_output_kmer(peptide_list, expr_lists, k)
    true_output = ['MTH\t1\t913.33', 'THA\t1\t870.0', 'HAW\t1\t740.0']
    assert c == true_output


def test_get_concat_peptide():
    front_coord = '10;19;25;33'
    back_coord = '27;36;44;53'
    front_peptide = ''
    back_peptide = 'MGF'
    strand = '+'
    concat_pep = get_concat_peptide(front_coord,back_coord,front_peptide,back_peptide,strand)
    assert concat_pep == ''

    front_peptide = 'MGF'
    back_peptide = ''
    strand = '+'
    concat_pep = get_concat_peptide(front_coord,back_coord,front_peptide,back_peptide,strand)
    assert concat_pep == ''

    front_peptide = 'EDM'
    back_peptide = 'DMF'
    strand = '+'
    concat_pep = get_concat_peptide(front_coord,back_coord,front_peptide,back_peptide,strand)
    assert concat_pep == 'EDMF'

    # neg case
    front_coord = '35;43;20;29'
    back_coord = '18;26;17;13'
    strand = '-'
    front_peptide = 'EDM'
    back_peptide = 'DMF'
    concat_pep = get_concat_peptide(front_coord,back_coord,front_peptide,back_peptide,strand)
    assert concat_pep == 'EDMF'

