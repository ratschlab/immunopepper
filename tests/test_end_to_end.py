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



@pytest.mark.parametrize("test_type, test_id,case,mutation_mode", [
    ['_base', '1', 'pos', 'ref'],
    ['_base', '1', 'pos', 'germline'],
    ['_base', '1', 'pos', 'somatic'],
    ['_base', '1', 'pos', 'somatic_and_germline'],
    ['_base', '1', 'neg', 'ref'],
    ['_base', '1', 'neg', 'germline'],
    ['_base', '1', 'neg', 'somatic'],
    ['_base', '1', 'neg', 'somatic_and_germline'],
    ['_insertion', '2', 'pos', 'ref'],
    ['_insertion', '2', 'pos', 'germline'],
    ['_insertion', '2', 'pos', 'somatic'],
    ['_insertion', '2', 'pos', 'somatic_and_germline'],
    ['_insertion', '2', 'neg', 'ref'],
    ['_insertion', '2', 'neg', 'germline'],
    ['_insertion', '2', 'neg', 'somatic'],
    ['_insertion', '2', 'neg', 'somatic_and_germline']
])


def test_end_to_end_ref(test_type, test_id, case, mutation_mode, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_type), 'data')
    sample_dir = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_type), 'test{}{}'.format(test_id,case))

    out_dir = str(tmpdir)

    my_args = ['--samples', 'test{}{}'.format(test_id,case),
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
               '--mutation_mode', mutation_mode]

    main_immuno.main(main_immuno.parse_arguments(my_args))

    _assert_files_equal(
        os.path.join(sample_dir, '{}_peptides_gt.fa'.format(mutation_mode)),
        os.path.join(out_dir, 'test{}{}'.format(test_id, case), '{}_peptides.fa'.format(mutation_mode)))

    _assert_files_equal(
        os.path.join(sample_dir, '{}_metadata_gt.tsv'.format(mutation_mode)),
        os.path.join(out_dir, 'test{}{}'.format(test_id, case), '{}_metadata.tsv.gz'.format(mutation_mode)))

