from immunopepper import immunopepper
import os

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


def test_end_to_end_build_mouse(tmpdir, mutation_mode):
    data_dir = os.path.join(os.path.dirname(__file__), 'mouse_usecase')
    out_dir = str(tmpdir)
    #sample_dir_build = os.path.join(os.path.dirname(__file__), 'test{}'.format(test_id),'diff','{}'.format(case),'test{}{}'.format(test_id,case))

    my_args_build = ['build',
               '--samples', 'ERR2130621',
               '--output-dir', out_dir,
               '--splice-path',os.path.join(data_dir,'ImmunoPepper_usecase.pickle'),
               '--count-path', os.path.join(data_dir,'ImmunoPepper_usecase.count.hdf5'),
               '--ann-path', os.path.join(data_dir,'ImmunoPepper_usecase.gtf'),
               '--ref-path', os.path.join(data_dir,'GRCm38.p6.genome.fa'),
               '--germline', os.path.join(data_dir,'ImmunoPepper_usecase.vcf'),
               '--somatic', os.path.join(data_dir,'ImmunoPepper_usecase.maf'),
                '--mutation-mode', mutation_mode,
                '--kmer', '9', '--parallel', '4']

    my_args = my_args_build
    immunopepper.split_mode(my_args)


### Mouse Test
tmpdir = '/Users/laurieprelot/Documents/Projects/tmp_kmer'
mutation_mode ='somatic_and_germline'
test_end_to_end_build_mouse(tmpdir, mutation_mode)

### Human Test
# test_id ='1'
# case = 'pos'
# mutation_mode ='germline'
# test_end_to_end_build(test_id, case, mutation_mode, tmpdir)