"""
Unit test for imunopepper/mutations.py
"""
import os

import numpy as np
import pytest

from spladder.classes.gene import Gene

import immunopepper.mutations as mutations
from immunopepper.namedtuples import Mutation
import immunopepper.ip as ip

test_data_dir = os.path.join(os.path.dirname(__file__), 'data')


class TestGetMutatedSequence:
    """ Unit-tests for :meth:`mutations.get_mutated_sequence` """

    def test_apply_mutations(self):
        ref_seq = 'GTAATGTGTAAGATGACGCACGCATGGTGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTTC'
        gt_mut_seq2 = 'GTAATGTGTAAGATGACGCACGCATA' + 'C' + 'TGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTTC'
        gt_mut_seq3 = 'GTAATGTGTAAGATGACGCACGCATA' + 'G' + 'TGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTTC'
        mut_dict = {}
        mut_dict[10] = {'mut_base': '*', 'ref_base': 'A'}
        mut_dict[25] = {'mut_base': 'A', 'ref_base': 'G'}
        mut_dict[26] = {'mut_base': 'C', 'ref_base': 'G'}
        mut_seq1 = mutations._apply_mutations(ref_seq[45:], 45, 46, mut_dict)  # test unclear mut_base
        assert mut_seq1 == ref_seq[45:]
        mut_seq2 = mutations._apply_mutations(ref_seq[25:], 25, 27, mut_dict)  # [25,27) include 26
        assert mut_seq2 == gt_mut_seq2[25:]
        mut_seq3 = mutations._apply_mutations(ref_seq[25:], 25, 26, mut_dict)  # [25,26) not include 26
        assert mut_seq3 == gt_mut_seq3[25:]

    def test_get_mutated_sequence(self):
        fasta_file = os.path.join(test_data_dir, 'small_genome.fa')

        ref_seq = 'GTAATGTGTAAGATGACGCACGCATGGTGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTTC'
        mutated_seq = 'GTAATGTGTAAGATGACGCACGCAT' + 'AC' + 'TGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTTC'
        mut_dict = {}
        mut_dict[9] = {'mut_base': 'G', 'ref_base': 'G'}  # outside range
        mut_dict[20] = {'mut_base': '*', 'ref_base': 'A'}
        mut_dict[35] = {'mut_base': 'A', 'ref_base': 'G'}
        mut_dict[36] = {'mut_base': 'C', 'ref_base': 'G'}
        mut_dict[75] = {'mut_base': 'T', 'ref_base': 'G'}  # outside range

        # the sub-sequence read from the fasta file is ref_seq
        seq_dict = mutations.get_mutated_sequence(fasta_file, 'chr1', 10, 74, mut_dict)  # test unclear mut_base
        assert ref_seq == seq_dict['ref']
        assert mutated_seq == seq_dict['background']

    def test_get_mutated_sequence_beginning_and_end(self):
        """ Test that mutations are correctly applied when they are right at the beginning and end """
        fasta_file = os.path.join(test_data_dir, 'small_genome.fa')

        ref_seq = 'GTAATGTGTAAGATGACGCACGCATGGTGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTTC'
        mutated_seq = 'A' + 'TAATGTGTAAGATGACGCACGCATGGTGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTT' + 'G'
        mut_dict = {}
        mut_dict[10] = {'mut_base': 'A', 'ref_base': 'G'}
        mut_dict[20] = {'mut_base': '*', 'ref_base': 'A'}
        mut_dict[73] = {'mut_base': 'G', 'ref_base': 'G'}

        # the sub-sequence read from the fasta file is ref_seq
        seq_dict = mutations.get_mutated_sequence(fasta_file, 'chr1', 10, 74, mut_dict)  # test unclear mut_base
        assert ref_seq == seq_dict['ref']
        assert mutated_seq == seq_dict['background']


class TestLoadMutations:
    @pytest.fixture
    def basic_args(self):
        return ['build',
                '--output-samples', 'sample1, sample2',
                '--splice-path', 'this_splicegraph',
                '--output-dir', 'this_output_dir',
                '--ann-path', 'this_ann_path',
                '--ref-path', 'this_ref_path']

    def test_reference(self, basic_args):
        """ Tests that when using --build-mode ref, the mutation dicts are empty """
        args = ip.parse_arguments(basic_args)
        mutation = mutations.load_mutations(args, args.output_samples)
        assert mutation.somatic_dict == {}
        assert mutation.germline_dict == {}

    def test_germline(self, basic_args):
        """ Tests that an error is thrown when requesting germline mutation mode, but --germline is missing. """
        my_args2 = basic_args + ['--germline', os.path.join(test_data_dir, 'test1pos.vcf')]
        args = ip.parse_arguments(my_args2)
        mutation = mutations.load_mutations(args, args.output_samples)
        assert mutation.somatic_dict == {}
        assert mutation.germline_dict == {
            ('test1pos', 'X'): {14: {'ref_base': 'G', 'mut_base': 'C', 'qual': '100', 'filter': 'PASS'}}}

    def test_somatic(self, basic_args):
        """ Tests that an error is thrown when requesting somatic mutation mode, but --somatic is missing. """
        my_args2 = basic_args + ['--somatic', os.path.join(test_data_dir, 'test1pos.maf')]
        args = ip.parse_arguments(my_args2)
        mutation = mutations.load_mutations(args, args.output_samples)
        assert mutation.germline_dict == {}
        assert mutation.somatic_dict == {('test1pos', 'X'): {
            38: {'ref_base': 'A', 'mut_base': 'G', 'strand': '+', 'variant_Classification': 'Silent',
                 'variant_Type': 'SNP'}}}

    def test_germline_and_somatic(self, basic_args):
        """ Tests that an error is thrown when requesting somatic mutation mode, but --somatic is missing. """
        my_args2 = basic_args + ['--germline', os.path.join(test_data_dir, 'test1pos.vcf'),
                                 '--somatic', os.path.join(test_data_dir, 'test1pos.maf')]
        args = ip.parse_arguments(my_args2)
        mutation = mutations.load_mutations(args, args.output_samples)
        assert mutation.somatic_dict == {('test1pos', 'X'): {
            38: {'ref_base': 'A', 'mut_base': 'G', 'strand': '+', 'variant_Classification': 'Silent',
                 'variant_Type': 'SNP'}}}
        assert mutation.germline_dict == {
            ('test1pos', 'X'): {14: {'ref_base': 'G', 'mut_base': 'C', 'qual': '100', 'filter': 'PASS'}}}


class TestGetSubMutations:
    @pytest.fixture
    def mutation(self):
        germline_dict = {('test1pos', '1'): {135: {'mut_base': 'G', 'ref_base': 'C'}},
                         ('test1neg', '1'): {136: {'mut_base': 'G', 'ref_base': 'C'}}}
        somatic_dict = {('test1pos', '1'): {28: {'mut_base': 'G', 'ref_base': 'C'}},
                        ('test2neg', '1'): {29: {'mut_base': 'G', 'ref_base': 'C'}}}
        mutation = Mutation(mode='somatic_and_germline', germline_dict=germline_dict, somatic_dict=somatic_dict)
        return mutation

    def test_get_sub_mutations(self, mutation):
        sample = 'test1pos'
        chrm = '1'

        sub_mutation = mutations.get_sub_mutations(mutation, sample, chrm)
        assert sub_mutation.somatic_dict == mutation.somatic_dict[('test1pos', '1')]
        assert sub_mutation.germline_dict == mutation.germline_dict[('test1pos', '1')]
        assert sub_mutation.mode == mutation.mode

    def test_inexistent_sample(self, mutation):
        # (test1neg,'1') not in the keys of maf_dict
        sample = 'inexistent'
        chrm = '1'
        sub_mutation = mutations.get_sub_mutations(mutation, sample, chrm)
        assert sub_mutation.somatic_dict == {}
        assert sub_mutation.germline_dict == {}
        assert sub_mutation.mode == mutation.mode

    def test_inexistent_chromosome(self, mutation):
        # (test1pos, 2) does exists neither in vcf nor maf
        sample = 'test1pos'
        chrm = '2'
        sub_mutation = mutations.get_sub_mutations(mutation, sample, chrm)
        assert sub_mutation.somatic_dict == {}
        assert sub_mutation.germline_dict == {}
        assert sub_mutation.mode == mutation.mode

    def test_somatic_only(self, mutation):
        """ Test a sample for which we only have somatic mutations """
        sample = 'test2neg'
        chrm = '1'

        sub_mutation = mutations.get_sub_mutations(mutation, sample, chrm)
        assert sub_mutation.somatic_dict == mutation.somatic_dict[('test2neg', '1')]
        assert sub_mutation.germline_dict == {}
        assert sub_mutation.mode == mutation.mode


def test_get_mut_comb():
    exon_to_mutations = {1: [100, 200], 2: [300]}
    exond_ids = [1, 2]
    assert [np.nan, (100,), (200,), (300,), (100, 200), (100, 300), (200, 300),
            (100, 200, 300)] == mutations.get_mut_comb(exon_to_mutations, exond_ids)


def test_exon_to_mutations():
    gene = Gene('ENSMUSG00000025902.13', 1000, 1200, 'chr1', '+')
    gene.splicegraph.vertices = np.array([[4491718, 4493771, 4495135], [4491720, 4493863, 4495155]], dtype='int')
    gene.splicegraph.edges = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype='int')

    assert {0: [4491719], 1: [], 2: [], np.nan: []} == mutations.exon_to_mutations(gene, [100, 200, 4491719])


class TestReadMafVcf:
    """ Tests reading and caching MAF and VCF files """

    def test_reading_vcf_h5(self):
        vcf_dict_default_heter_code0 = mutations.parse_mutation_from_vcf('somatic',
                                                                         os.path.join(test_data_dir, 'test1vcf.h5'),
                                                                         target_samples=['test1pos', 'test1neg'])
        vcf_dict_heter_code2 = mutations.parse_mutation_from_vcf('somatic', os.path.join(test_data_dir, 'test1vcf.h5'),
                                                                 target_samples=['test1pos', 'test1neg'], heter_code=2)
        assert len(vcf_dict_default_heter_code0) == 2
        assert vcf_dict_default_heter_code0['test1neg', 'X'][135] == {'mut_base': 'G', 'ref_base': 'C'}
        assert vcf_dict_default_heter_code0['test1pos', 'X'][14] == {'mut_base': 'C', 'ref_base': 'G'}
        assert vcf_dict_heter_code2['test1pos', 'X'][135] == {'mut_base': 'G', 'ref_base': 'C'}
        assert vcf_dict_heter_code2['test1neg', 'X'][14] == {'mut_base': 'C', 'ref_base': 'G'}
        assert vcf_dict_heter_code2['test1neg', 'X'][135] == {'mut_base': 'G', 'ref_base': 'C'}

    def test_cache_maf(self):
        pass
