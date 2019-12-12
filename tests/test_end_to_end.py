import gzip
import os

import pytest


data_dir = os.path.join(os.path.dirname(__file__), 'data')
groundtruth_dir = os.path.join(os.path.dirname(__file__), 'test1')

from immunopepper import main_immuno


def _assert_files_equal(expected_path, actual_path):
    def o(f):
        if f.endswith('.gz'):
            return gzip.open(f, 'rb')
        else:
            return open(f, 'r')

    with o(expected_path) as e:
        with o(actual_path) as a:
            afile = a.read()
            efile = e.read()
            assert afile == efile

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

# test build mode
def test_end_to_end_build(test_id, case, mutation_mode, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_id), 'data')
    out_dir = str(tmpdir)
    sample_dir_build = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_id),'build','{}'.format(case),'test{}{}'.format(test_id,case))

    my_args_build = ['build','--samples', 'test{}{}'.format(test_id,case),
               '--output-dir', out_dir,
               '--splice-path',
               '{}/{}graph/spladder/genes_graph_conf3.merge_graphs.pickle'.format(
                   data_dir, case),
               '--count-path',
               '{}/{}graph/spladder/genes_graph_conf3.merge_graphs.count.hdf5'.format(
                   data_dir, case),
               '--ann-path', '{}/test{}{}.gtf'.format(data_dir, test_id, case),
               '--ref-path', '{}/test{}{}.fa'.format(data_dir, test_id, case),
               '--germline', '{}/test{}{}.vcf'.format(data_dir, test_id, case),
               '--somatic', '{}/test{}{}.maf'.format(data_dir, test_id, case),
               '--mutation-mode', mutation_mode,
               '--kmer', '4']

    my_args = my_args_build
    sample_dir = sample_dir_build
    main_immuno.split_mode(my_args)
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
    #kmer
    _assert_files_equal(
        os.path.join(sample_dir, '{}_back_kmer.txt'.format(mutation_mode)),
        os.path.join(out_dir, 'test{}{}'.format(test_id, case), '{}_back_kmer.txt'.format(mutation_mode)))
    _assert_files_equal(
        os.path.join(sample_dir, '{}_junction_kmer.txt'.format(mutation_mode)),
        os.path.join(out_dir, 'test{}{}'.format(test_id, case), '{}_junction_kmer.txt'.format(mutation_mode)))
    _assert_files_equal(
        os.path.join(sample_dir, 'gene_expression_detail.tsv'),
        os.path.join(out_dir, 'test{}{}'.format(test_id, case), 'gene_expression_detail.tsv'))

# test
@pytest.mark.parametrize("test_id,case", [
    ['1', 'pos'],
    ['1', 'neg']
])
def test_end_to_end_makebg(test_id, case,tmpdir):
    sample_dir = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_id),'make_bg','{}'.format(case))
    out_dir = str(tmpdir)
    output_file_name = os.path.join(out_dir,'{}_integrated_background_kmer.txt'.format(case))
    back_file_dir = os.path.join(os.path.dirname(__file__),'test{}'.format(test_id),'build',case,'test{}{}'.format(test_id,case))
    bg_file_list = [os.path.join(back_file_dir,'{}_back_kmer.txt'.format(mode)) for mode in ['ref','somatic','somatic_and_germline','germline']]
    my_args_makebg = ['make_bg','--kmer-files']+bg_file_list+['--output-file-path', output_file_name]+['--output-dir', out_dir]

    main_immuno.split_mode(my_args_makebg)
    _assert_files_equal(
        os.path.join(sample_dir, '{}_integrated_background_kmer.txt'.format(case)),
        os.path.join(out_dir,'{}_integrated_background_kmer.txt'.format(case)))

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
def test_end_to_end_diff(test_id, case, tmpdir,mutation_mode):
    groundtruth_diff_dir = os.path.join(groundtruth_dir,'diff','{}'.format(case))
    out_dir = str(tmpdir)
    junction_kmer_file_path = os.path.join(groundtruth_dir, 'build',
                                           case, 'test{}{}'.format(test_id, case),
                                           '{}_junction_kmer.txt'.format(mutation_mode))

    bg_kmer_file_path = os.path.join(groundtruth_dir, 'make_bg',
                                     '{}'.format(case), '{}_integrated_background_kmer.txt'.format(case))

    output_file_path = os.path.join(out_dir,'{}_{}_junction_kmer_with_bg.txt'.format(case, mutation_mode))

    my_args_diff = ['diff', '--junction-kmer-file', junction_kmer_file_path,
                   '--bg-file-path', bg_kmer_file_path,
                   '--output-file-path', output_file_path,
                   '--output-dir', out_dir]

    main_immuno.split_mode(my_args_diff)
    _assert_files_equal(
        os.path.join(groundtruth_diff_dir, '{}_{}_junction_kmer_with_bg.txt'.format(case,mutation_mode)),
        os.path.join(out_dir,'{}_{}_junction_kmer_with_bg.txt'.format(case,mutation_mode)))

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
def test_end_to_end_filter(test_id, case, tmpdir,mutation_mode):
    groundtruth_filter_dir = os.path.join(groundtruth_dir,'filter','{}'.format(case))
    out_dir = str(tmpdir)
    junction_kmer_file_path = os.path.join(groundtruth_dir, 'diff',
                                           case,'{}_{}_junction_kmer_with_bg.txt'.format(case, mutation_mode))

    # cross_junction filter
    output_filtered_file_name = 'cj_{}_{}_junction_kmer_with_bg.txt'.format(case, mutation_mode)
    output_file_path = os.path.join(out_dir,output_filtered_file_name)
    my_args = ['filter', '--junction-kmer-tsv-path', junction_kmer_file_path,
               '--output-file-path', output_file_path,
               '--output-dir', out_dir,
               '--cross-junction']
    main_immuno.split_mode(my_args)
    _assert_files_equal(
        os.path.join(groundtruth_filter_dir, output_filtered_file_name),
        os.path.join(out_dir,output_filtered_file_name))

    # junction expression 0 filter
    output_filtered_file_name = 'junc_expr_0_{}_{}_junction_kmer_with_bg.txt'.format(case, mutation_mode)
    output_file_path = os.path.join(out_dir,output_filtered_file_name)
    my_args = ['filter', '--junction-kmer-tsv-path', junction_kmer_file_path,
               '--output-file-path', output_file_path,
               '--output-dir', out_dir,
               '--junc-expr']
    main_immuno.split_mode(my_args)
    _assert_files_equal(
        os.path.join(groundtruth_filter_dir, output_filtered_file_name),
        os.path.join(out_dir, output_filtered_file_name))

    # segment expression 1300 filter
    output_filtered_file_name = 'seg_expr_1300_{}_{}_junction_kmer_with_bg.txt'.format(case, mutation_mode)
    output_file_path = os.path.join(out_dir,output_filtered_file_name)
    my_args = ['filter', '--junction-kmer-tsv-path', junction_kmer_file_path,
               '--output-file-path', output_file_path,
               '--output-dir', out_dir,
               '--seg-expr',
               '--seg-expr-thresh','1300']
    main_immuno.split_mode(my_args)
    _assert_files_equal(
        os.path.join(groundtruth_filter_dir, output_filtered_file_name),
        os.path.join(out_dir, output_filtered_file_name))

    # segment expression 1300 filter and junction expression 500
    output_filtered_file_name = 'junc_expr_500_seg_expr_1300_cj_{}_{}_junction_kmer_with_bg.txt'.format(case, mutation_mode)
    output_file_path = os.path.join(out_dir,output_filtered_file_name)
    my_args = ['filter', '--junction-kmer-tsv-path', junction_kmer_file_path,
               '--output-file-path', output_file_path,
               '--output-dir', out_dir,
               '--seg-expr',
               '--seg-expr-thresh', '1300',
               '--junc-expr',
               '--junc-expr-thresh', '500']
    main_immuno.split_mode(my_args)
    _assert_files_equal(
        os.path.join(groundtruth_filter_dir, output_filtered_file_name),
        os.path.join(out_dir, output_filtered_file_name))
