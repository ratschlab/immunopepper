import gzip
import os

import pytest


data_dir = os.path.join(os.path.dirname(__file__), 'data')
sample_dir = os.path.join(os.path.dirname(__file__), 'test1')
from scripts import main_immuno


def _assert_files_equal(expected_path, actual_path):
    def o(f):
        if f.endswith('.gz'):
            return gzip.open(f, 'rb')
        else:
            return open(f, 'r')

    with o(expected_path) as e:
        with o(actual_path) as a:
            assert e.read() == a.read()



@pytest.mark.parametrize("test_type, case, mutation_mode", [
    ['_base', 'pos', 'ref'],
    ['_base', 'pos', 'germline'],
    ['_base', 'pos', 'somatic'],
    ['_base', 'pos', 'somatic_and_germline'],
    ['_base', 'neg', 'ref'],
    ['_base', 'neg', 'germline'],
    ['_base', 'neg', 'somatic'],
    ['_base', 'neg', 'somatic_and_germline'],
    ['_insertion', 'pos', 'ref'],
    ['_insertion', 'pos', 'germline'],
    ['_insertion', 'pos', 'somatic'],
    ['_insertion', 'pos', 'somatic_and_germline'],
    ['_insertion', 'neg', 'ref'],
    ['_insertion', 'neg', 'germline'],
    ['_insertion', 'neg', 'somatic'],
    ['_insertion', 'neg', 'somatic_and_germline']
])


def test_end_to_end_ref(test_type, case, mutation_mode, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_type), 'data')
    sample_dir = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_type), 'test{}{}'.format(test_type,case))

    out_dir = str(tmpdir)

    my_args = ['--samples', 'test{}{}'.format(test_type,case),
               '--output_dir', out_dir,
               '--splice_path',
               '{}/{}graph/spladder/genes_graph_conf3.merge_graphs.pickle'.format(
                   data_dir, case),
               '--count_path',
               '{}/{}graph/spladder/genes_graph_conf3.merge_graphs.count.hdf5'.format(
                   data_dir, case),
               '--ann_path', '{}/test{}{}.gtf'.format(data_dir, test_type, case),
               '--ref_path', '{}/test{}{}.fa'.format(data_dir, test_type, case),
               '--vcf_path', '{}/test{}{}.vcf'.format(data_dir, test_type, case),
               '--maf_path', '{}/test{}{}.maf'.format(data_dir, test_type, case),
               '--mutation_mode', mutation_mode,
               '--is_filter']

    main_immuno.main(main_immuno.parse_arguments(my_args))

    _assert_files_equal(
        os.path.join(sample_dir, '{}_peptides_gt.fa'.format(mutation_mode)),
        os.path.join(out_dir, 'test{}{}'.format(test_type, case), '{}_peptides.fa'.format(mutation_mode)))

    _assert_files_equal(
        os.path.join(sample_dir, '{}_metadata_gt.tsv'.format(mutation_mode)),
        os.path.join(out_dir, 'test{}{}'.format(test_type, case), '{}_metadata.tsv.gz'.format(mutation_mode)))

