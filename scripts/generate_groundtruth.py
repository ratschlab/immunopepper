import os
from immunopepper import main_immuno
test_id = 1
case_list = ['pos','neg']
mutation_list = ['ref','somatic','germline','somatic_and_germline']

test_file_dir = '/Users/jiayu/PycharmProjects/CBM_RA/projects2018_immunopepper/tests'

for case in case_list:
    for mutation_mode in mutation_list:
        data_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'data')
        out_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'groundtruth_kmer',
                                       '{}'.format(case))

        my_args = ['--samples', 'test{}{}'.format(test_id, case),
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
        main_immuno.main(main_immuno.parse_arguments(my_args))

# generate filtered junction output
for case in case_list:
    for mutation_mode in mutation_list:
        data_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'data')
        out_dir = os.path.join(test_file_dir, 'test{}'.format(test_id), 'groundtruth_junction',
                                       '{}'.format(case))

        my_args = ['--samples', 'test{}{}'.format(test_id, case),
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
                        '--filter_redundant']
        main_immuno.main(main_immuno.parse_arguments(my_args))

