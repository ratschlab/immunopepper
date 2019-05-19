####
# This script is used to generate the neo_kmer for immunopepper whose version is before 'write_pep_on_the_fly' branch.
# In the original version, the junction kmer and concat kmer are written in different files.
# also the junction_count is not included in kmer file (only exist in metadata file).
# This script will be deprecated once new version of immunopepper is used
import os
import pandas as pd
import numpy as np
from functools import reduce
from immunopepper.constant import NOT_EXIST

def get_unique_kmer(file,zero_filter=True,cross_junction=True):
    kmer_df = pd.read_csv(file,header=None,names=['kmer','gene_name','expr','flag'],sep='\t')
    kmer_df = kmer_df.dropna()
    kmer_df = kmer_df.astype({'expr':'float16'})
    kmer_df['flag'] = kmer_df['flag'].astype('bool')
    if cross_junction:
        kmer_df = kmer_df[kmer_df['flag']]
    if zero_filter:
        kmer_df = kmer_df[kmer_df['expr']>0]
    return kmer_df

def modify_gene_name(concat_kmer_df,gene_id_dict):
    def name_transfer(old_name):
        id = int(old_name.split('_')[0][4:])
        gene_name = gene_id_dict[id]
        other_info = old_name.split('_')[1:]
        new_name_list = [gene_name] + other_info
        new_name = '_'.join(new_name_list)
        return new_name
    concat_kmer_df['gene_name'] = concat_kmer_df['gene_name'].apply(lambda x: name_transfer(x))
    return concat_kmer_df

def add_expression_count(kmer_df,meta_df):
    def get_expression_count(gene_name_str,junction_series):
        items = gene_name_str.split('_')
        gene_name = items[0]
        if len(items) == 4:
            junction1 = items[1]+','+items[2]
            count1 = junction_series[(gene_name,junction1)]
            junction2 = items[2]+','+items[3]
            count2 = junction_series[(gene_name,junction2)]
            if count1 != NOT_EXIST and count2 != NOT_EXIST:
                count = np.mean((float(count1),float(count2)))
            else:
                count = 0
        else:
            assert len(items) == 3
            junction1 = items[1]+','+items[2]
            count1 = junction_series[(gene_name,junction1)]
            if count1 != NOT_EXIST:
                count = float(count1)
            else:
                count = 0
        return count
    junction_series = meta_df.groupby(['gene_name', 'vertex_idx'])['junction_expr'].agg('max') # index by junction
    kmer_df['junction_expr'] = kmer_df['gene_name'].apply(lambda x: get_expression_count(x,junction_series))
    return kmer_df

def integrate_concate_kmer(data_dir,gene_id_dict):
    mode_prefix = ['germline','somatic','somatic_and_germline']
    kmer_df_list = []
    for mode in mode_prefix:
        concat_file_name = os.path.join(data_dir,mode+'_concat_kmer.txt')
        junction_file_name = os.path.join(data_dir,mode+'_junction_kmer.txt')
        meta_file_name = os.path.join(data_dir,mode+'_metadata.tsv.gz')
        meta_df = pd.read_csv(meta_file_name, sep='\t')
        concat_kmer_df = get_unique_kmer(concat_file_name,zero_filter=False,cross_junction=False)

        # remove false positive cross-junction concat_kmer
        # those not appear in the junction_kmer_df are all cross-junction kmer
        junction_kmer_df_all = get_unique_kmer(junction_file_name,zero_filter=False,cross_junction=False)
        concat_cross_junction_kmer_set = set(concat_kmer_df['kmer']).difference(junction_kmer_df_all['kmer'])
        concat_cross_junction_kmer_df = concat_kmer_df.loc[concat_kmer_df['kmer'].isin(concat_cross_junction_kmer_set)]
        concat_cross_junction_kmer_df = modify_gene_name(concat_cross_junction_kmer_df,gene_id_dict)

        junction_kmer_df = get_unique_kmer(junction_file_name,zero_filter=False,cross_junction=True)
        final_junction_kmer_df = pd.concat((concat_cross_junction_kmer_df,junction_kmer_df))

        final_junction_kmer_df = final_junction_kmer_df.groupby('kmer')[['expr', 'gene_name']].agg('max')
        final_junction_kmer_df['mode'] = mode
        final_junction_kmer_df = add_expression_count(final_junction_kmer_df,meta_df)

        # only take kmer with junction expression count greater than 0
        final_junction_kmer_df = final_junction_kmer_df[final_junction_kmer_df['junction_expr']>0]
        kmer_df_list.append(final_junction_kmer_df)
    final_junction_kmer_df_all_mode = reduce(lambda x,y: pd.concat((x,y)), kmer_df_list)
    annotaion_background_kmer_list = []
    for mode in ['germline','somatic']:
        back_file_name = os.path.join(data_dir,mode+'_back_kmer.txt')
        junction_kmer_df_all = get_unique_kmer(back_file_name, zero_filter=False, cross_junction=False)
        annotaion_background_kmer_list.extend(list(junction_kmer_df_all['kmer']))

    # subtract background kmer
    final_junction_kmer_df_all_mode_subtract_background = final_junction_kmer_df_all_mode[np.logical_not(final_junction_kmer_df_all_mode.index.isin(annotaion_background_kmer_list))]
    return final_junction_kmer_df_all_mode_subtract_background



gene_id_dict = np.load('gene_id_dict.npy').item()
result_dir = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun'
all_dir = list(filter(lambda dir_name: dir_name.startswith('TCGA'),os.listdir(result_dir)))

dir_name = 'TCGA-24-1430'
mode = 'germline'
data_dir = os.path.join(result_dir, dir_name)

file_num_list = [len(os.listdir(os.path.join(result_dir,dir_name))) for dir_name in all_dir]
for dir_name in all_dir:
    print(dir_name)
    data_dir = os.path.join(result_dir,dir_name)
    final_kmer_df = integrate_concate_kmer(data_dir,gene_id_dict)
    result_file = os.path.join(data_dir,'integrate_kmer')
    final_kmer_df.to_csv(result_file,sep='\t')
