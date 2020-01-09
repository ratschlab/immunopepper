import os
from immunopepper import main_immuno
test_id = 1
case_list = ['pos','neg']
mutation_list = ['ref','somatic','germline','somatic_and_germline']

test_file_dir = '../tests'

# build mode
for case in case_list:
    for mutation_mode in mutation_list:
        data_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'data')
        out_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'build',case)

        my_args = ['build','--samples', 'test{}{}'.format(test_id, case),
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
        main_immuno.split_mode(my_args)

# make_bg mode
for case in case_list:
    output_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'make_bg',case)
    output_file_name = os.path.join(out_dir,'{}_integrated_background_kmer.txt'.format(case))
    back_file_dir = os.path.join(test_file_dir,'test{}'.format(test_id),'build',case,'test{}{}'.format(test_id,case))
    bg_file_list = [os.path.join(back_file_dir,'{}_back_kmer.txt'.format(mode)) for mode in ['ref','somatic','somatic_and_germline','germline']]
    my_args = ['make_bg','--kmer-files']+bg_file_list+['--output-file-path',output_file_name]+['--output-dir',output_dir]
    main_immuno.split_mode(my_args)

# diff mode
for case in case_list:
    for mutation_mode in mutation_list:
        junction_kmer_file_path = os.path.join(test_file_dir, 'test{}'.format(test_id), 'build',
                                               case, 'test{}{}'.format(test_id, case),'{}_junction_kmer.txt'.format(mutation_mode))

        bg_kmer_file_path = os.path.join(test_file_dir, 'test{}'.format(test_id), 'make_bg',
                                         '{}'.format(case), '{}_integrated_background_kmer.txt'.format(case))
        output_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'diff','{}'.format(case))
        output_file_path = os.path.join(output_dir, '{}_{}_junction_kmer_with_bg.txt'.format(case,mutation_mode))
        output_file_path_remove_bg = os.path.join(output_dir, '{}_{}_junction_kmer_without_bg.txt'.format(case,mutation_mode))

        my_args = ['diff', '--junction-kmer-file', junction_kmer_file_path,
                   '--bg-file-path', bg_kmer_file_path,
                   '--output-file-path', output_file_path,
                   '--output-dir',output_dir]
        main_immuno.split_mode(my_args)

        my_args = ['diff', '--junction-kmer-file', junction_kmer_file_path,
                   '--bg-file-path', bg_kmer_file_path,
                   '--output-file-path', output_file_path_remove_bg,
                   '--remove-bg',
                   '--output-dir', output_dir
                   ]
        main_immuno.split_mode(my_args)
# filter mode
for case in case_list:
    for mutation_mode in mutation_list:
        output_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'filter',
                                        '{}'.format(case))
        junction_kmer_tsv_path = os.path.join(test_file_dir, 'test{}'.format(test_id), 'diff',
                                               '{}'.format(case), '{}_{}_junction_kmer_with_bg.txt'.format(case,mutation_mode))
        meta_file_path = os.path.join(test_file_dir, 'test{}'.format(test_id), 'build',
                                               '{}'.format(case), 'test{}{}'.format(test_id,case),'{}_metadata.tsv.gz'.format(mutation_mode))
        junction_kmer_tsv_file_name = '{}_{}_junction_kmer_with_bg.txt'.format(case,mutation_mode)
        my_args = ['filter', '--junction-kmer-tsv-path', junction_kmer_tsv_path,
                   '--output-file-path', os.path.join(output_dir,'cj_'+junction_kmer_tsv_file_name),
                   '--cross-junction',
                   '--output-dir', output_dir]
        main_immuno.split_mode(my_args)

        my_args = ['filter', '--junction-kmer-tsv-path', junction_kmer_tsv_path,
                   '--output-file-path', os.path.join(output_dir,'junc_expr_0_'+junction_kmer_tsv_file_name),
                   '--junc-expr',
                   '--output-dir', output_dir]
        main_immuno.split_mode(my_args)

        my_args = ['filter', '--junction-kmer-tsv-path', junction_kmer_tsv_path,
                   '--output-file-path', os.path.join(output_dir,'seg_expr_1300_'+junction_kmer_tsv_file_name),
                   '--seg-expr',
                   '--seg-expr-thresh',str(1300),
                   '--output-dir',output_dir]
        main_immuno.split_mode(my_args)

        my_args = ['filter', '--junction-kmer-tsv-path', junction_kmer_tsv_path,
                   '--output-file-path', os.path.join(output_dir,'seg_expr_1300_peptide_annotated_1_junction_annotated_1_stop_1_isolated_1_'+junction_kmer_tsv_file_name),
                   '--meta-file-path', meta_file_path,
                   '--seg-expr',
                   '--seg-expr-thresh',str(1300),
                   '--peptide-annotated',str(1),
                   '--junction-annotated',str(1),
                   '--has-stop-codon',str(1),
                   '--is-isolated',str(0),
                   '--output-dir',output_dir]
        main_immuno.split_mode(my_args)


        my_args = ['filter', '--junction-kmer-tsv-path', junction_kmer_tsv_path,
                   '--output-file-path', os.path.join(output_dir,'junc_expr_500_seg_expr_1300_cj_'+junction_kmer_tsv_file_name),
                   '--cross-junction',
                   '--junc-expr',
                   '--junc-expr-thresh',str(500),
                   '--seg-expr',
                   '--seg-expr-thresh', str(1300),
                   '--output-dir',output_dir]
        main_immuno.split_mode(my_args)

