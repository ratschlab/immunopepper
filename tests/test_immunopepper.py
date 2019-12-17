import gzip
import os
import pickle

import Bio.SeqIO as BioIO
import pytest
import numpy as np

from immunopepper.immuno_mutation import apply_germline_mutation,construct_mut_seq_with_str_concat,get_mutation_mode_from_parser,Mutation,get_sub_mutation_tuple
from immunopepper.immuno_preprocess import preprocess_ann, genes_preprocess, \
    parse_mutation_from_vcf, parse_mutation_from_maf
from immunopepper.utils import get_sub_mut_dna,get_concat_peptide,convert_namedtuple_to_str,check_chr_consistence,get_idx,create_libsize
from immunopepper.io_utils import load_pickled_graph
from immunopepper.main_immuno import parse_arguments
from immunopepper.immuno_model import create_output_kmer
from immunopepper.immuno_nametuple import Coord,OutputBackground,OutputKmer
from immunopepper.constant import NOT_EXIST
from immunopepper.immuno_filter import get_junction_anno_flag
data_dir = os.path.join(os.path.dirname(__file__), 'test1','data')


@pytest.fixture
def load_gene_data():
    f = open(os.path.join(data_dir, 'posgraph','spladder',
                          'genes_graph_conf3.merge_graphs.pickle'), 'rb')
    ann_path = os.path.join(data_dir, 'test1pos.gtf')
    ref_path = os.path.join(data_dir, 'test1pos.fa')

    (graph_data, graph_meta) = load_pickled_graph(f)  # cPickle.load(f)
    genetable,chr_set = preprocess_ann(ann_path)
    interesting_chr = list(map(str, range(1, 23))) + ["X", "Y", "MT"]
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
    mutation_dic_vcf = parse_mutation_from_vcf(vcf_path=vcf_path,h5_sample_list=['test1pos'])
    mutation_dic_maf = parse_mutation_from_maf(maf_path=maf_path)

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
                                          mutation_sub_dict=mutation_sub_dic_vcf)
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
        coord = Coord(vlist[0], vlist[1], vlist[2], vlist[3])
        sub_dna = get_sub_mut_dna(ref_seq, coord, variant_comb[i],
                                  mutation_sub_dic_maf, strand[i])
        assert sub_dna == groundtruth[i]


def test_reading_gtf_and_gff3_file():
    gff3_path = os.path.join(data_dir,'small.gencode.v29.gff3')
    gtf_path = os.path.join(data_dir, 'small.gencode.v29.gtf')
    gene_table_gtf,_ = preprocess_ann(gtf_path)
    gene_table_gff,_ = preprocess_ann(gff3_path)
    assert gene_table_gff.gene_to_ts == gene_table_gtf.gene_to_ts
    assert gene_table_gtf.ts_to_cds.keys() == gene_table_gff.ts_to_cds.keys()
    assert gene_table_gtf.ts_to_cds['ENST00000335137.4'] == [(69090, 70005, 0)]
    assert gene_table_gff.ts_to_cds['ENST00000335137.4'] == [(69090, 70008, 0)] # include stop codon in gff3


def test_reading_gtf_and_gff_file():
    gff_path = os.path.join(data_dir,'small.gencode.v19.gff')
    gtf_path = os.path.join(data_dir, 'small.gencode.v19.gtf')
    gene_table_gtf,_ = preprocess_ann(gtf_path)
    gene_table_gff,_ = preprocess_ann(gff_path)
    assert gene_table_gff.gene_to_ts == gene_table_gtf.gene_to_ts
    assert gene_table_gtf.ts_to_cds.keys() == gene_table_gff.ts_to_cds.keys()
    assert gene_table_gff.ts_to_cds['ENST00000335137.3'] == [(69090, 70006, 0)]
    assert gene_table_gtf.ts_to_cds['ENST00000335137.3'] == [(69090, 70005, 0)]


def test_reading_vcf_h5():
    vcf_dict_default_heter_code0 = parse_mutation_from_vcf(os.path.join(data_dir,'test1vcf.h5'),h5_sample_list=['test1pos','test1neg'])
    vcf_dict_heter_code2 = parse_mutation_from_vcf(os.path.join(data_dir,'test1vcf.h5'),h5_sample_list=['test1pos','test1neg'],heter_code=2)
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
    basic_args = ['build',
                  '--samples','this_sample',
                  '--splice-path','this_splicegraph',
                  '--output-dir','this_output_dir',
                  '--ann-path','this_ann_path',
                  '--ref-path','this_ref_path']
    my_args1 = basic_args+[
               '--germline', os.path.join(data_dir,'test1pos.vcf'),
               '--somatic', os.path.join(data_dir,'test1pos.maf'),
               '--mutation-mode', 'somantic']  # bad mutation mode
    args = parse_arguments(my_args1)
    try:
        get_mutation_mode_from_parser(args)
    except SystemExit:
        assert 1
    my_args2 = basic_args+['--germline', os.path.join(data_dir,'test1pos.vcf'),
               '--mutation-mode', 'somatic']  # mismatch mutation mode and input files
    args = parse_arguments(my_args2)
    try:
        get_mutation_mode_from_parser(args)
    except SystemExit:
        assert 1

def test_get_sub_mutation_tuple():
    germline_dict = {('test1pos','1'):{135:{'mut_base': 'G', 'ref_base': 'C'}},('test1neg','1'):{136:{'mut_base': 'G', 'ref_base': 'C'}}}
    somatic_dict = {('test1pos','1'):{28:{'mut_base': 'G', 'ref_base': 'C'}},('test2neg','1'):{29:{'mut_base': 'G', 'ref_base': 'C'}}}
    sample = 'test1pos'
    chrm = '1'
    mutation = Mutation(mode='somatic_and_germline',germline_mutation_dict=germline_dict,somatic_mutation_dict=somatic_dict)

    sub_mutation = get_sub_mutation_tuple(mutation,sample,chrm)
    assert sub_mutation.somatic_mutation_dict == mutation.somatic_mutation_dict[('test1pos','1')]
    assert sub_mutation.germline_mutation_dict == mutation.germline_mutation_dict[('test1pos','1')]
    assert sub_mutation.mode == mutation.mode

    # (test1neg,'1') not in the keys of maf_dict
    sample = 'test1neg'
    try:
        sub_mutation = get_sub_mutation_tuple(mutation, sample, chrm)
    except SystemExit:
        assert 1

    # (test1pos, 2) does exists neither in vcf nor maf
    sample = 'test1pos'
    chrm = '2'
    try:
        sub_mutation = get_sub_mutation_tuple(mutation, sample, chrm)
    except SystemExit:
        assert 1

def test_create_output_kmer():
    k = 3
    peptide = OutputBackground('1','MTHAW')
    expr_lists = [(8,1000),(1,220),(6,0)] # test 0 expression
    c = create_output_kmer(peptide, expr_lists, k)
    true_output = [OutputKmer('MTH','1',913.33,False,NOT_EXIST), OutputKmer('THA','1',580.0,False,NOT_EXIST), OutputKmer('HAW','1',246.67,False,NOT_EXIST)]
    assert c == true_output
    expr_lists = [(8,1000),(1,220),(0,0)] # test 0 expression
    c = create_output_kmer(peptide, expr_lists, k)
    true_output = [OutputKmer('MTH','1',913.33,False,NOT_EXIST), OutputKmer('THA','1',870.0,False,NOT_EXIST), OutputKmer('HAW','1',740.0,False,NOT_EXIST)]
    assert c == true_output


def test_get_concat_peptide():
    front_coord = Coord(10,19,25,33)
    back_coord = Coord(27,36,44,53)
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
    front_coord = Coord(35,43,20,29)
    back_coord = Coord(18,26,17,13)
    strand = '-'
    front_peptide = 'EDM'
    back_peptide = 'DMF'
    concat_pep = get_concat_peptide(front_coord,back_coord,front_peptide,back_peptide,strand)
    assert concat_pep == 'EDMF'

def test_convert_namedtuple_to_str():
    other_pep_field_list = ['id', 'new_line', 'peptide']
    back_pep1 = OutputBackground(1,'EDMHG')
    back_pep2 = OutputBackground(2,'')
    back_pep3 = OutputBackground(3,'KKQ')
    back_pep_list = [back_pep1,back_pep2,back_pep3]
    result = [convert_namedtuple_to_str(back_pep,other_pep_field_list)+'\n' for back_pep in back_pep_list]
    expected_result = ['1\nEDMHG\n', '2\n\n', '3\nKKQ\n']
    assert result == expected_result

    other_pep_field_list = ['kmer','id','expr','is_cross_junction','junction_count']
    kmer_pep1 = OutputKmer('','GENE0_1_2',NOT_EXIST,False,NOT_EXIST)
    kmer_pep2 = OutputKmer('AQEB','GENE0_1_3',20,True,25)
    kmer_pep_list = [kmer_pep1,kmer_pep2]
    result = [convert_namedtuple_to_str(kmer_pep,other_pep_field_list)+'\n' for kmer_pep in kmer_pep_list]
    expected_result = ['\tGENE0_1_2\t.\tFalse\t.\n', 'AQEB\tGENE0_1_3\t20\tTrue\t25\n']
    assert result == expected_result

def test_get_junction_ann_flag():
    junction_flag = np.zeros((4,4))
    junction_flag[1,2] = 1
    junction_flag[2,3] = 1
    vertex_id_tuple = (1,2,3)
    assert get_junction_anno_flag(junction_flag,vertex_id_tuple) == [1,1]
    vertex_id_tuple = (0,2,3)
    assert get_junction_anno_flag(junction_flag,vertex_id_tuple) == [0,1]
    vertex_id_tuple = (0,1,3)
    assert get_junction_anno_flag(junction_flag,vertex_id_tuple) == [0,0]
    vertex_id_tuple = (1,2)
    assert get_junction_anno_flag(junction_flag,vertex_id_tuple) == 1

def test_check_chr_consistence(load_gene_data, load_mutation_data):
    graph_data, ref_seq, gene_cds_begin_dict = load_gene_data
    mutation_dic_vcf, mutation_dic_maf = load_mutation_data
    mutation = Mutation(None,mutation_dic_maf,mutation_dic_vcf) # vcf_dict[('test1pos','X')]
    ann_path = os.path.join(data_dir, 'test1pos.gtf')
    genetable,chr_set = preprocess_ann(ann_path)

    check_chr_consistence(chr_set,mutation,graph_data)

    # conflict between chr_set and vcf_dict
    try:
        chr_set = set(['chrX'])
        check_chr_consistence(chr_set, mutation, graph_data)
    except SystemExit:
        assert 1

    # conflict between chr_set and gene.chr
    try:
        graph_data[0].chr = 'chrX'
        check_chr_consistence(chr_set, mutation, graph_data)
    except SystemExit:
        assert 1


def test_get_idx():
    sample_idx_table = {}
    sample_idx_table['test1pos'] = 100
    sample_idx_table['test1neg'] = 101
    sample = 'test1pos'
    gene_idx = 0
    Idx = get_idx(sample_idx_table,sample,gene_idx)
    assert Idx.gene == gene_idx
    assert Idx.sample == sample_idx_table['test1pos']

    # the inquiry sample name not in sample_idx_table
    Idx = get_idx(sample_idx_table, 'test2', gene_idx)
    assert Idx.gene == gene_idx
    assert Idx.sample == None

    # the sample_idx_table does not exist
    Idx = get_idx(None, 'test2', gene_idx)
    assert Idx.gene == gene_idx
    assert Idx.sample == None

def test_create_libsize():
    fp = 'temp.txt'
    expr_distr_dict = {'test1pos':['.'],'test1neg':[2,10,11]}
    libsize_count_dict = create_libsize(expr_distr_dict,fp,debug=True)
    assert libsize_count_dict == {'test1neg':(10.5,23)}
