import gzip
import os

import pytest


data_dir = os.path.join(os.path.dirname(__file__), 'data')
sample_dir = os.path.join(os.path.dirname(__file__), 'test1')
from immunopepper import main_immuno


def _assert_files_equal(expected_path, actual_path):
    def o(f):
        if f.endswith('.gz'):
            return gzip.open(f, 'rb')
        else:
            return open(f, 'r')

    with o(expected_path) as e:
        with o(actual_path) as a:
            assert e.read() == a.read()

@pytest.mark.parametrize("test_id,case,mutation_mode", [
    ['1', 'pos', 'ref'],
    ['1', 'pos', 'germline'],
    ['1', 'pos', 'somatic'],
    ['1', 'pos', 'somatic_and_germline'],
    ['1', 'neg', 'ref'],
    ['1', 'neg', 'germline'],
    ['1', 'neg', 'somatic'],
    ['1', 'neg', 'somatic_and_germline']
])


def test_end_to_end_ref(test_id, case, mutation_mode, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_id), 'data')
    out_dir = str(tmpdir)
    sample_dir_kmer = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_id),'groundtruth_kmer','{}'.format(case),'test{}{}'.format(test_id,case))
    sample_dir_junction = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_id),'groundtruth_junction','{}'.format(case),'test{}{}'.format(test_id,case))

    my_args_kmer = ['--samples', 'test{}{}'.format(test_id,case),
               '--output_dir', out_dir,
               '--splice_path',
               '{}/{}graph/spladder/genes_graph_conf3.merge_graphs.pickle'.format(
                   data_dir, case),
               '--count_path',
               '{}/{}graph/spladder/genes_graph_conf3.merge_graphs.count.hdf5'.format(
                   data_dir, case),
               '--ann_path', '{}/test{}{}.gtf'.format(data_dir, test_id, case),
               '--ref_path', '{}/test{}{}.fa'.format(data_dir, test_id, case),
               '--vcf_path', '{}/test{}{}.vcf'.format(data_dir, test_id, case),
               '--maf_path', '{}/test{}{}.maf'.format(data_dir, test_id, case),
               '--mutation_mode', mutation_mode,
               '--output_silence',
                '--kmer', '4']
    my_args_junction = my_args_kmer[:-2]+['--filter_redundant']

    for my_args,sample_dir in zip([my_args_kmer],[sample_dir_kmer]):
        main_immuno.main(main_immuno.parse_arguments(my_args))
        _assert_files_equal(
            os.path.join(sample_dir, '{}_metadata.tsv.gz'.format(mutation_mode)),
            os.path.join(out_dir, 'test{}{}'.format(test_id, case), '{}_metadata.tsv.gz'.format(mutation_mode)))
        # peptide
        _assert_files_equal(
            os.path.join(sample_dir, '{}_peptides.fa'.format(mutation_mode)),
            os.path.join(out_dir, 'test{}{}'.format(test_id, case), '{}_peptides.fa'.format(mutation_mode)))
        _assert_files_equal(
            os.path.join(sample_dir, '{}_back_peptides.fa'.format(mutation_mode)),
            os.path.join(out_dir, 'test{}{}'.format(test_id, case), '{}_back_peptides.fa'.format(mutation_mode)))
        _assert_files_equal(
            os.path.join(sample_dir, '{}_concat_peptides.fa'.format(mutation_mode)),
            os.path.join(out_dir, 'test{}{}'.format(test_id, case), '{}_concat_peptides.fa'.format(mutation_mode)))
        #kmer
        _assert_files_equal(
            os.path.join(sample_dir, '{}_back_kmer.txt'.format(mutation_mode)),
            os.path.join(out_dir, 'test{}{}'.format(test_id, case), '{}_back_kmer.txt'.format(mutation_mode)))
        _assert_files_equal(
            os.path.join(sample_dir, '{}_concat_kmer.txt'.format(mutation_mode)),
            os.path.join(out_dir, 'test{}{}'.format(test_id, case), '{}_concat_kmer.txt'.format(mutation_mode)))
        _assert_files_equal(
            os.path.join(sample_dir, '{}_junction_kmer.txt'.format(mutation_mode)),
            os.path.join(out_dir, 'test{}{}'.format(test_id, case), '{}_junction_kmer.txt'.format(mutation_mode)))
