import gzip
import os

import pytest

from immunopepper import main_immuno

data_dir = os.path.join(os.path.dirname(__file__), 'data')
sample_dir = os.path.join(os.path.dirname(__file__), 'test1')


def _assert_files_equal(expected_path, actual_path):
    def o(f):
        if f.endswith('.gz'):
            return gzip.open(f, 'rb')
        else:
            return open(f, 'r')

    with o(expected_path) as e:
        with o(actual_path) as a:
            assert e.read() == a.read()


@pytest.mark.parametrize("vcf_path,maf_path,prefix,", [
    [None, None, 'ref'],
    ['test1.vcf', None, 'germline'],
    [None, 'test1.maf', 'somatic'],
    ['test1.vcf', 'test1.maf', 'somatic_and_germline']
])
def test_end_to_end_ref(vcf_path, maf_path, prefix, tmpdir):
    out_dir = str(tmpdir)

    my_args = ['--samples', 'test1',
               '--output_dir', out_dir,
               '--splice_path',
               '{}/spladder/genes_graph_conf3.merge_graphs.pickle'.format(
                   data_dir),
               '--count_path',
               '{}/spladder/genes_graph_conf3.merge_graphs.count_0ts.hdf5'.format(
                   data_dir),
               '--ann_path', '{}/test1.gtf'.format(data_dir),
               '--ref_path', '{}/test1.fa'.format(data_dir)]

    if vcf_path:
        my_args.extend(['--vcf_path', os.path.join(data_dir, vcf_path)])

    if maf_path:
        my_args.extend(['--maf_path', os.path.join(data_dir, maf_path)])

    main_immuno.main(main_immuno.parse_arguments(my_args))

    _assert_files_equal(
        os.path.join(sample_dir, '{}_peptides_gt.fa'.format(prefix)),
        os.path.join(out_dir, 'test1', '{}_peptides.fa'.format(prefix)))

    _assert_files_equal(
        os.path.join(sample_dir, '{}_metadata_gt.tsv'.format(prefix)),
        os.path.join(out_dir, 'test1', '{}_metadata.tsv.gz'.format(prefix)))
