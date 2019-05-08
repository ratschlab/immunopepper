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
        junction_kmer_df = get_unique_kmer(junction_file_name,junction_zero_filter=True,cross_junction=True)
        final_junction_kmer_df = junction_kmer_df.groupby('kmer')[['seg_expr', 'gene_name','junction_expr']].agg('max')
        final_junction_kmer_df['mode'] = mode
        kmer_df_list.append(final_junction_kmer_df)
    final_junction_kmer_df_all_mode = reduce(lambda x,y: pd.concat((x,y)), kmer_df_list)
    annotation_background_kmer_list = []
    for mode in ['germline','somatic']:
        back_file_name = os.path.join(data_dir,mode+'_back_kmer.txt')
        junction_kmer_df_all = get_unique_kmer(back_file_name, junction_zero_filter=False, cross_junction=False)
        annotation_background_kmer_list.extend(list(junction_kmer_df_all['kmer']))

    # subtract background kmer
    final_junction_kmer_df_all_mode_subtract_background = final_junction_kmer_df_all_mode[np.logical_not(final_junction_kmer_df_all_mode.index.isin(annotaion_background_kmer_list))]
    return final_junction_kmer_df_all_mode_subtract_background



gene_id_dict = np.load('gene_id_dict.npy').item()
result_dir = '/cluster/work/grlab/projects/TCGA/immunopepper_tcga_rerun_fly'
all_dir = list(filter(lambda dir_name: dir_name.startswith('TCGA'),os.listdir(result_dir)))

dir_name = 'TCGA-24-1430'
mode = 'germline'
data_dir = os.path.join(result_dir, dir_name)

file_num_list = [len(os.listdir(os.path.join(result_dir,dir_name))) for dir_name in all_dir]
for dir_name in all_dir:
    print(dir_name)
    data_dir = os.path.join(result_dir,dir_name)
    final_kmer_df = integrate_concate_kmer(data_dir)
    result_file = os.path.join(data_dir,'integrate_kmer')
    final_kmer_df.to_csv(result_file,sep='\t')
