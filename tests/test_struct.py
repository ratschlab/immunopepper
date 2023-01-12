from immunopepper import ip
import os


def test_end_to_end_build_mouse(tmpdir, mutation_mode, is_parallel=True, graph_cross_sample=False):
    data_dir = os.path.join(os.path.dirname(__file__), 'mouse_usecase')
    out_dir = str(tmpdir)
    #sample_dir_build = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_id),'diff','{}'.format(case),'test{}{}'.format(test_id,case))
    my_args_build = ['build',
                #'--kmer-database', os.path.join(data_dir, 'uniprot_test.csv'),
                #'--libsize-extract',
               #'--disable-process-libsize',
               '--output-samples', 'ERR2130621', #ERR2130621',# 'ENCSR000BZG'
               '--mutation-sample', 'ERR2130621',
               #'--pickle-samples', 'ERR2130621',
               '--output-dir', out_dir,
               '--splice-path', os.path.join(data_dir,'ImmunoPepper_usecase.pickle'),
               '--count-path', os.path.join(data_dir,'ImmunoPepper_usecase.count.hdf5'),
               '--ann-path', os.path.join(data_dir,'ImmunoPepper_usecase.gtf'),
               '--ref-path', os.path.join(data_dir,'GRCm38.p6.genome.fa'),
               '--kmer', '34',
               '--batch-size', '1',
               '--output-fasta',
               #'--all-read-frames',
               #'--process-num', '1',
               '--verbose', '1',
               #'--min-pep-len', '2',
               #'--process-chr', 'chr2',
               #'--genes-interest', '/Users/laurieprelot/Documents/Projects/tmp_kmer/restrict_genes_test/genes_of_interest.tsv'
                  ]
    if mutation_mode == 'germline' or mutation_mode == 'somatic_and_germline':
        my_args_build.extend(['--germline', os.path.join(data_dir,'ImmunoPepper_usecase.vcf')])
    if mutation_mode == 'somatic' or mutation_mode == 'somatic_and_germline':
        my_args_build.extend(['--somatic', os.path.join(data_dir,'ImmunoPepper_usecase.maf')])
    if is_parallel:
        my_args_build.extend(['--parallel', '4'])
    if graph_cross_sample:
        my_args_build.extend(['--cross-graph-exp' ])
    my_args = my_args_build
    ip.split_mode(my_args)


def test_end_to_end_samplespecif(sample, tmpdir, kmer, mutation_mode):
    out_dir = os.path.join(tmpdir, sample)
    bg_kmer_file_path = os.path.join(out_dir, 'integrated_background_kmer.pq.gz')
    junction_kmer_file_path = os.path.join(out_dir,
                                           '{}_sample_{}mer.pq.gz'.format(mutation_mode, kmer))
    output_file_path = os.path.join(out_dir,
                                    '{}_sample_{}mer_with_bg.pq.gz'.format(mutation_mode, kmer))

    bg_file_list = [os.path.join(out_dir, '{}_annot_{}mer.pq.gz'.format(mode, kmer)) for mode in
                    ['ref', 'somatic', 'somatic_and_germline', 'germline']]

    my_args= ['samplespecif', '--kmer-files'] + bg_file_list +  [
        '--output-dir', out_dir, '--compressed',
        '--junction-kmer-file', junction_kmer_file_path,
        '--bg-file-path', bg_kmer_file_path,
        '--output-file-path', output_file_path,
        '--remove-bg']

    ip.split_mode(my_args)


def test_end_to_end_cancerspecif_mx():

    basedir = "/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test"
    my_args =["cancerspecif",
              "--cores", "2",
              "--mem-per-core", "6000",
              "--kmer", "9",
              #"--path-cancer-libsize",os.path.join(basedir,'cancer_no_ct_var', 'libsize_cancer.tsv'),
              #"--path-normal-libsize", os.path.join(basedir, 'normal', 'libsize_normals_top20'),
              #"--path-cancer-matrix-segm", os.path.join(basedir, 'cancer', 'ref_graph_kmer_SegmExpr_top20_n20_overlap.pq.gz'),
              "--path-cancer-matrix-edge", os.path.join(basedir, 'cancer', 'simple_foreground_flag.gz'), #TODO moce to proper test folder
              "--path-normal-matrix-segm", os.path.join(basedir, 'normal',  'simple_background_flag_bool.gz'),  os.path.join(basedir, 'normal', 'nested' ,'simple_background_flag_copy_flag_bool.gz'),
              "--path-normal-matrix-edge", os.path.join(basedir, 'normal', 'nested', 'simple_background_overlay1_flag_bool.gz'),
              # "--path-normal-matrix-edge", os.path.join(basedir, 'normal', 'ref_graph_kmer_SegmExpr_top20_n20_overlap.pq.gz'),
              #"--path-normal-kmer-list", os.path.join(basedir, 'normal', 'simple_annotation.pq'),
              # "--path-normal-kmer-list", os.path.join(basedir, 'normal', 'simple_background_precompute.pq'),os.path.join(basedir, 'normal', 'simple_annotation.pq'),
              #"--path-normal-kmer-list", os.path.join(basedir, 'normal', 'simple_background_precompute.tsv'),os.path.join(basedir, 'normal', 'simple_background_precompute_bis.tsv'),
              '--ids-cancer-samples',  "TCGA-13-1497-01A-01" , "TCGA-24-1103-01A-01",
              "--output-dir", os.path.join(basedir, 'filter_out'),
              "--output-count", os.path.join(basedir, 'filter_out', "collect.txt"),
              '--expr-high-limit-normal', "0.0",
              '--cohort-expr-support-normal', "0",
              "--n-samples-lim-normal", "0",
              '--sample-expr-support-cancer', "4",
              '--cohort-expr-support-cancer', "20",
              "--n-samples-lim-cancer", "1",
              "--filterNeojuncCoord", "C",
              #"--filterAnnotatedRF", "",
              # "--tot-batches", "4",
              # "--batch-id", "0",
              "--tag-prefix", 'G_',
              "--interm-dir-norm", "/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test/test_inter",
              "--interm-dir-canc", "/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test/test_inter_c",
              #"--expression-fields-c", 'segment_expr', 'junction_expr',
              #"--tissue-grp-files", "/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test/normal/tissue_grps/dummy_BRCA.txt",
              #'/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test/normal/tissue_grps/dummy_OV.txt',
              #"--whitelist-normal", "/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test/normal/tissue_grps/dummy_BRCA.txt" ,
              '--whitelist-cancer', '/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test/cancer/whitelist_cancer.txt',
              "--uniprot", "/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test/uniprot",
              #"--path-normal-kmer-list", "/Users/laurieprelot/Documents/Projects/tmp_kmer/filter_test/filter_out/normals_merge-segm-edge_max_uniq_expr-in-3-samples-with-0.1-normalized-cts.tsv",
              "--parallelism", "3",
              "--out-partitions", "2",
              "--mut-cancer-samples", "ref", "ref"]
              #"--statistical"]
    print(my_args)
    ip.split_mode(my_args)





### Mouse Test
tmpdir = '/Users/laurieprelot/Documents/Projects/tmp_kmer'
mutation_mode ='germline'
#pr = cProfile.Profile()
#pr.enable()
#test_end_to_end_build_mouse(tmpdir, mutation_mode, is_parallel=False) #TODO add back
#test_end_to_end_samplespecif('ERR2130621', tmpdir, "9", mutation_mode) # TEST DEPRECATED
#for mutation_mode in ['ref', 'germline', 'somatic', 'somatic_and_germline']:
#     test_end_to_end_makebg('ERR2130621', tmpdir, "9")
#     test_end_to_end_diff(tmpdir, 'ERR2130621', "9", mutation_mode)
test_end_to_end_cancerspecif_mx()
#pr.disable()
#pr.dump_stats(os.path.join(tmpdir, 'cProfile.pstats'))


