import os
from immunopepper import main_immuno
test_id = 1
case_list = ['pos','neg']
mutation_list = ['ref','somatic','germline','somatic_and_germline']

test_file_dir = '/Users/jiayu/PycharmProjects/CBM_RA/projects2018_immunopepper/tests'

# build mode
for case in case_list:
    for mutation_mode in mutation_list:
        data_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'data')
        out_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'build',case)

        my_args = ['build','--samples', 'test{}{}'.format(test_id, case),
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
                        '--kmer', '4']
        main_immuno.split_mode(my_args)

# make_bg mode
for case in case_list:
    out_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'make_bg',case)
    output_file_name = os.path.join(out_dir,'{}_integrated_background_kmer.txt'.format(case))
    back_file_dir = os.path.join(test_file_dir,'test{}'.format(test_id),'build',case,'test{}{}'.format(test_id,case))
    bg_file_list = [os.path.join(back_file_dir,'{}_back_kmer.txt'.format(mode)) for mode in ['ref','somatic','somatic_and_germline','germline']]
    my_args = ['make_bg','--kmer_files_list']+bg_file_list+['--output_file_path',output_file_name]
    main_immuno.split_mode(my_args)

# diff mode
for case in case_list:
    for mutation_mode in mutation_list:
        junction_kmer_file_path = os.path.join(test_file_dir, 'test{}'.format(test_id), 'build',
                                               case, 'test{}{}'.format(test_id, case),'{}_junction_kmer.txt'.format(mutation_mode))

        bg_kmer_file_path = os.path.join(test_file_dir, 'test{}'.format(test_id), 'make_bg',
                                         '{}'.format(case), '{}_integrated_background_kmer.txt'.format(case))

        output_file_path = os.path.join(test_file_dir, 'test{}'.format(test_id), 'diff',
                                        '{}'.format(case), '{}_{}_junction_kmer_with_bg.txt'.format(case,mutation_mode))

        output_file_path_remove_bg = os.path.join(test_file_dir, 'test{}'.format(test_id), 'diff',
                                        '{}'.format(case), '{}_{}_junction_kmer_without_bg.txt'.format(case,mutation_mode))

        my_args = ['diff', '--junction_kmer_file', junction_kmer_file_path,
                   '--bg_file_path', bg_kmer_file_path,
                   '--output_file_path', output_file_path]
        main_immuno.split_mode(my_args)

        my_args = ['diff', '--junction_kmer_file', junction_kmer_file_path,
                   '--bg_file_path', bg_kmer_file_path,
                   '--output_file_path', output_file_path_remove_bg,
                   '--remove_bg']
        main_immuno.split_mode(my_args)
# filter mode
for case in case_list:
    for mutation_mode in mutation_list:
        output_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'filter',
                                        '{}'.format(case))
        junction_kmer_tsv_path = os.path.join(test_file_dir, 'test{}'.format(test_id), 'diff',
                                               '{}'.format(case), '{}_{}_junction_kmer_with_bg.txt'.format(case,mutation_mode))

        junction_kmer_tsv_file_name = '{}_{}_junction_kmer_with_bg.txt'.format(case,mutation_mode)
        my_args = ['filter', '--junction_kmer_tsv_path', junction_kmer_tsv_path,
                   '--output_file_path', os.path.join(output_dir,'cj_'+junction_kmer_tsv_file_name),
                   '--cross_junction']
        main_immuno.split_mode(my_args)

        my_args = ['filter', '--junction_kmer_tsv_path', junction_kmer_tsv_path,
                   '--output_file_path', os.path.join(output_dir,'junc_expr_0_'+junction_kmer_tsv_file_name),
                   '--junc_expr']
        main_immuno.split_mode(my_args)

        my_args = ['filter', '--junction_kmer_tsv_path', junction_kmer_tsv_path,
                   '--output_file_path', os.path.join(output_dir,'seg_expr_1300_'+junction_kmer_tsv_file_name),
                   '--seg_expr',
                   '--seg_expr_thre',str(1300)]
        main_immuno.split_mode(my_args)

        my_args = ['filter', '--junction_kmer_tsv_path', junction_kmer_tsv_path,
                   '--output_file_path', os.path.join(output_dir,'junc_expr_500_seg_expr_1300_cj_'+junction_kmer_tsv_file_name),
                   '--cross_junction',
                   '--junc_expr',
                   '--junc_expr_thre',str(500),
                   '--seg_expr',
                   '--seg_expr_thre', str(1300)]
        main_immuno.split_mode(my_args)

