import pstats
import os
import profile

from immunopepper import immunopepper


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
                '--kmer', '9']

    my_args = my_args_build
    immunopepper.split_mode(my_args)

def main():
    ### Mouse Test
    tmpdir = '/Users/laurieprelot/Documents/Projects/tmp_kmer'
    mutation_mode ='somatic'
    test_end_to_end_build_mouse(tmpdir, mutation_mode)


profile.run('main()', '~/my_profile_tmp')

stats = pstats.Stats('~/my_profile_tmp')
os.remove('~/my_profile_tmp')
stats.sort_stats('cumtime').print_stats('immunopepper')