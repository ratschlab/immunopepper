from immunopepper import immunopepper
import os
import cProfile
import pathlib


def test_end_to_end_build(test_id, case, mutation_mode, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_id), 'data')
    out_dir = str(tmpdir)
    sample_dir_build = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_id),'diff','{}'.format(case),'test{}{}'.format(test_id,case))

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
                '--gtex-junction-path', '{}/{}graph/spladder/genes_graph_conf3.test1{}.junction.hdf5'.format(data_dir, case, case),
                '--mutation-mode', mutation_mode,
                '--kmer', '4']

    my_args = my_args_build
    sample_dir = sample_dir_build
    immunopepper.split_mode(my_args)


def test_end_to_end_build_mouse(tmpdir, mutation_mode, is_parallel=True):
    data_dir = os.path.join(os.path.dirname(__file__), 'mouse_usecase')
    out_dir = str(tmpdir)
    #sample_dir_build = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_id),'diff','{}'.format(case),'test{}{}'.format(test_id,case))
    my_args_build = ['build',
               '--samples', 'ENCSR000BZG', 'ERR2130621',#'ERR2130621','ENCSR000BZG'
               '--output-dir', out_dir,
               '--splice-path',os.path.join(data_dir,'ImmunoPepper_usecase.pickle'),
               '--count-path', os.path.join(data_dir,'ImmunoPepper_usecase.count.hdf5'),
               '--ann-path', os.path.join(data_dir,'ImmunoPepper_usecase.gtf'),
               '--ref-path', os.path.join(data_dir,'GRCm38.p6.genome.fa'),
               '--germline', os.path.join(data_dir,'ImmunoPepper_usecase.vcf'),
               '--somatic', os.path.join(data_dir,'ImmunoPepper_usecase.maf'),
                '--mutation-mode', mutation_mode,
                '--kmer', '9',
                '--batch-size', '1']
    if is_parallel:
        my_args_build.extend(['--parallel', '4'])

    my_args = my_args_build
    immunopepper.split_mode(my_args)

def test_end_to_end_makebg(sample, tmpdir, kmer):
    out_dir = os.path.join(tmpdir, sample)
    output_file_name = os.path.join(tmpdir, sample, 'integrated_background_kmer.pq.gz')
    back_file_dir = os.path.join(tmpdir, sample)
    bg_file_list = [os.path.join(back_file_dir, '{}_back_{}mer.pq.gz'.format(mode, kmer)) for mode in
                    ['ref', 'somatic', 'somatic_and_germline', 'germline']]
    my_args_makebg = ['make_bg', '--kmer-files'] + bg_file_list + ['--output-file-path', output_file_name] + [
        '--output-dir', out_dir, '--compressed']

    immunopepper.split_mode(my_args_makebg)

def test_end_to_end_diff(tmpdir, sample, kmer_length, mutation_mode):
    out_dir = str(tmpdir)
    junction_kmer_file_path = os.path.join(out_dir, sample, '{}_junction_{}mer.pq.gz'.format(mutation_mode, kmer_length))
    bg_kmer_file_path = os.path.join(out_dir, sample, 'integrated_background_kmer.pq.gz')
    output_file_path = os.path.join(out_dir, sample,'{}_junction_{}mer_with_bg.pq.gz'.format(mutation_mode, kmer_length))
    my_args_diff = ['diff', '--junction-kmer-file', junction_kmer_file_path,
                   '--bg-file-path', bg_kmer_file_path,
                   '--output-file-path', output_file_path,
                   '--output-dir', out_dir,
                   '--remove-bg']
    immunopepper.split_mode(my_args_diff)


def test_end_to_end_filter(tmpdir, sample, kmer_length, mutation_mode):
    out_dir = str(tmpdir)
    junction_kmer_file_path = os.path.join(out_dir, sample,
                            '{}_junction_{}mer_with_bg.pq.gz'.format(mutation_mode, kmer_length))
    meta_file_path = os.path.join(out_dir, sample,
                                  '{}_metadata.tsv.gz.pq'.format(mutation_mode))

    # segment expression 1300 filter and junction expression 500
    output_filtered_file_name = 'junc_expr_500_seg_expr_1300_cj_{}_junction_{}mer_with_bg.pq.gz'.format(mutation_mode, kmer_length)
    output_file_path = os.path.join(out_dir, sample, output_filtered_file_name)

    my_args = ['filter', '--junction-kmer-tsv-path', junction_kmer_file_path,
               '--output-file-path', output_file_path,
               '--output-dir', out_dir,
               '--seg-expr',
               '--seg-expr-thresh', '0.12',
               #'--junc-expr',
               #'--junc-expr-thresh', '0',
               '--meta-file-path', meta_file_path,
               '--peptide-annotated', '0',
               '--junction-annotated', '0',
               '--has-stop-codon', '0',
               '--is-in-junction-list', '1',
               '--is-isolated', '0']
    immunopepper.split_mode(my_args)

def test_end_to_end_crosscohort(tmpdir):
    my_args =["crosscohort",
              "--cores", "4",
              "--mem-per-core", "5000",
              "--mutation-modes","ref",
              "--kmer", "9",
              "--samples", "ERR2130621", "ENCSR000BZG",
              "--input-dir", tmpdir,
              "--output-dir", tmpdir,
              "--output-suffix", "test",
              "--compressed_inputs",
              "--skip-filegrouping"]
              #"--segment-count"]
    immunopepper.split_mode(my_args)

def mini_crosscohort():

    cancer_dir = "/Users/laurieprelot/Documents/Projects/tmp_kmer/dev_samples/cancer"
    my_args =["crosscohort",
              "--cores", "4",
              "--mem-per-core", "5000",
              "--mutation-modes","ref", "germline", "somatic", "somatic_and_germline",
              "--kmer", "9",
              "--samples", "TCGA-AO-A12D-01A-11", "TCGA-AR-A0TT-01A-31",
              "--input-dir", cancer_dir,
              "--output-dir", cancer_dir,
              "--output-suffix", "_test",
              "--compressed_inputs",
              "--remove-bg",
              "--skip-filegrouping"]
    immunopepper.split_mode(my_args)

def test_end_to_end_cancerspecif():

    basedir = "/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test"
    my_args =["cancerspecif",
              "--cores", "2",
              "--mem-per-core", "6000",
              "--kmer", "9",
              "--path-cancer-libsize",os.path.join(basedir,'cancer', 'expression_counts.libsize.tsv'),
              "--path-normal-libsize", os.path.join(basedir, 'normal', 'libsize_normals_top20'),
              "--paths-cancer-samples",
              "/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test/cancer_no_ct_var/TCGA-13-1497-01A-01/tmp_out_somatic_10000/somatic_junction_9mer_n20.pq.gz", "/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test/cancer_no_ct_var/TCGA-24-1103-01A-01/tmp_out_somatic_1000/somatic_junction_9mer_n20.pq.gz",
              "--path-normal-matrix-segm", os.path.join(basedir, 'normal', 'ref_graph_kmer_SegmExpr_top20_n20_overlap.pq.gz'),
              "--path-normal-matrix-edge", os.path.join(basedir, 'normal' 'ref_graph_kmer_JuncExpr.pq.gz'),
              "--output-dir", os.path.join(basedir, 'filter_out'),
              '--expr-high-limit-normal', "2.0",
              '--expr-limit-normal', "2.0",
              "--expr-n-limit", "1",
              "--statistical"]
    immunopepper.split_mode(my_args)

### Mouse Test
tmpdir = '/Users/laurieprelot/Documents/Projects/tmp_kmer'
mutation_mode ='ref'
#pr = cProfile.Profile()
#pr.enable()
#for mutation_mode in ['ref', 'germline', 'somatic', 'somatic_and_germline']:
#    test_end_to_end_build_mouse(tmpdir, mutation_mode, is_parallel=True) #TODO add back
#     test_end_to_end_makebg('ERR2130621', tmpdir, "9")
#     test_end_to_end_diff(tmpdir, 'ERR2130621', "9", mutation_mode)
#    test_end_to_end_filter(tmpdir, 'ERR2130621', "9", mutation_mode)
#test_end_to_end_crosscohort(tmpdir) #TODO add back
#mini_crosscohort()
test_end_to_end_cancerspecif()
#pr.disable()
#pr.dump_stats(os.path.join(tmpdir, 'cProfile.pstats'))


