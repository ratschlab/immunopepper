####
# This script is used for generating the neokmer (for CancerCell Paper)
# In the foreground part (line 49- line 64), we need only take the junction kmer and non-zero junction count kmer.
# After appling all this filters, we get a 'integrate_kmer' file in the corresponding sample subdirectory.
# The background part is easier to process since we take all the kmers so we just need a sort | uniq bash script.
# After getting the foreground kmer and background kmer, we will add a new column to the foreground kmer file
#   to indicate if the kmer is neo-kmer (not present in the background kmer). This part is in line 68 - line 95
# Once to be noticed is that, there are some kmer result for big graphs in `immunopepper_gtex_rerun_big_graph`,
#   and we also need to subtract them.
# This script is an updated version of `filter_kmer.py` which is compatible with the latest version of immunopepper.
import os
import pandas as pd
import numpy as np
from functools import reduce


def get_unique_kmer(file,file_type='background',junction_zero_filter=True,cross_junction=True,seg_zero_filter=False):
    if file_type == 'background':
        kmer_df = pd.read_csv(file, header=None, names=['kmer', 'gene_name', 'seg_expr', 'flag'],sep='\t')
    elif file_type == 'junction':
        kmer_df = pd.read_csv(file, header=None, names=['kmer', 'gene_name', 'seg_expr', 'flag', 'junction_expr'],sep='\t')
    else:
        print('file type {} is not recognized'.format(file_type))
        return 0
    kmer_df = kmer_df.dropna()
    kmer_df = kmer_df.astype({'seg_expr':'float16'})
    if cross_junction:
        kmer_df = kmer_df[kmer_df['flag']]
    if file_type == 'junction' and junction_zero_filter:
        kmer_df = kmer_df.astype({'junction_expr':'float16'})
        kmer_df = kmer_df[kmer_df['junction_expr']>0]
    if seg_zero_filter:
        kmer_df = kmer_df[kmer_df['seg_expr']>0]
    return kmer_df

def integrate_concate_kmer(data_dir):
    mode_prefix = ['germline','somatic','somatic_and_germline']
    kmer_df_list = []
    for mode in mode_prefix:
        junction_file_name = os.path.join(data_dir,mode+'_junction_kmer.txt')
        junction_kmer_df = get_unique_kmer(junction_file_name,file_type='junction',junction_zero_filter=True,cross_junction=True)
        final_junction_kmer_df = junction_kmer_df.groupby('kmer')[['seg_expr', 'gene_name','junction_expr']].agg('max')
        final_junction_kmer_df['mode'] = mode
        kmer_df_list.append(final_junction_kmer_df)
    final_junction_kmer_df_all_mode = reduce(lambda x,y: pd.concat((x,y)), kmer_df_list)
    annotation_background_kmer_list = []
    for mode in ['germline','somatic']:
        back_file_name = os.path.join(data_dir,mode+'_back_kmer.txt')
        junction_kmer_df_all = get_unique_kmer(back_file_name, file_type='background',junction_zero_filter=False, cross_junction=False)
        annotation_background_kmer_list.extend(list(junction_kmer_df_all['kmer']))

    # subtract background kmer
    final_junction_kmer_df_all_mode_subtract_background = final_junction_kmer_df_all_mode[np.logical_not(final_junction_kmer_df_all_mode.index.isin(annotation_background_kmer_list))]
    return final_junction_kmer_df_all_mode_subtract_background

result_dir = '/cluster/work/grlab/projects/TCGA/immunopepper_tcga_rerun_fly'
all_dir = list(filter(lambda dir_name: dir_name.startswith('TCGA'),os.listdir(result_dir)))

dir_name = 'TCGA-24-1430'
mode = 'germline'
data_dir = os.path.join(result_dir, dir_name)
for dir_name in all_dir:
    print(dir_name)
    data_dir = os.path.join(result_dir,dir_name)
    final_kmer_df = integrate_concate_kmer(data_dir)
    result_file = os.path.join(data_dir,'integrate_kmer')
    final_kmer_df.to_csv(result_file,sep='\t')


####
foreground_dir = '/cluster/work/grlab/projects/TCGA/immunopepper_tcga_rerun_fly'
background_dir = '/cluster/work/grlab/projects/TCGA/immunopepper_gtex_rerun'
background_big_dir = '/cluster/work/grlab/projects/TCGA/immunopepper_gtex_rerun_big_graph'

all_dir = list(filter(lambda dir_name: dir_name.startswith('TCGA'),os.listdir(background_dir)))
background_ref_file = os.path.join(background_dir, 'uniq_final_background_ref_kmer.txt')
background_ref_list = open(background_ref_file, 'r').read().strip().split('\n')

dir_name = 'TCGA-24-1430'
for dir_name in all_dir:
    print(dir_name)
    foreground_sample_dir = os.path.join(foreground_dir,dir_name)
    foreground_ref_file = os.path.join(foreground_sample_dir,'integrate_ref_kmer.txt')
    new_foreground_ref_file = os.path.join(foreground_sample_dir,'final_integrate_ref_kmer.txt')
    foreground_ref_kmer_df = pd.read_csv(foreground_ref_file,sep='\t')
    uniq_ref_kmer_set = set(foreground_ref_kmer_df['kmer']).difference(set(background_ref_list))
    foreground_ref_kmer_df['in_neo_flag'] = foreground_ref_kmer_df['kmer'].apply(lambda x:x in uniq_ref_kmer_set)
    foreground_ref_kmer_df.to_csv(new_foreground_ref_file,sep='\t')

    background_sample_dir = os.path.join(background_dir,dir_name)
    foreground_mut_file = os.path.join(foreground_sample_dir,'integrate_kmer')
    foreground_mut_kmer_df = pd.read_csv(foreground_mut_file,sep='\t')
    background_mut_file = os.path.join(background_sample_dir,'integrate_kmer.txt')
    new_foreground_mut_file = os.path.join(foreground_sample_dir,'final_integrate_mut_kmer.txt')
    background_mut_list = open(background_mut_file,'r').read().strip().split('\n')
    uniq_mut_kmer_set = set(foreground_mut_kmer_df['kmer']).difference(set(background_mut_list))
    foreground_mut_kmer_df['in_neo_flag'] = foreground_mut_kmer_df['kmer'].apply(lambda x:x in uniq_mut_kmer_set)
    foreground_mut_kmer_df.to_csv(new_foreground_mut_file,sep='\t')
