import os
import pandas as pd
import numpy as np
import sys
from immunopepper.immuno_print import print_memory_diags

def get_unique_kmer(file,zero_filter=True,cross_junction=True):
    kmer_df = pd.read_csv(file,header=None,names=['kmer','gene_name','expr','flag'],sep='\t')
    kmer_df = kmer_df.dropna()
    kmer_df = kmer_df.astype({'expr':'float16'})
    if cross_junction:
        kmer_df = kmer_df[np.logical_or(kmer_df['flag'] == 'True', kmer_df['flag'] == True)]
    if zero_filter:
        kmer_df = kmer_df[kmer_df['expr']>0]
    return kmer_df

##### part 2
# Gtex_variant
# get unique kmer set

def generate_unique_kmer_and_write(data_dir,mode='variant'):
    if mode=='variant':
        mode_prefix = ['germline','somatic','somatic_and_germline']
    elif mode == 'ref':
        mode_prefix = ['ref']
    else:
        print('mode can only choose between variant and ref')
    kmer_df_set = set()
    for mode in mode_prefix:
        concat_file_name = os.path.join(data_dir,mode+'_concat_kmer.txt')
        junction_file_name = os.path.join(data_dir,mode+'_junction_kmer.txt')
        concat_kmer_df = get_unique_kmer(concat_file_name,zero_filter=False,cross_junction=False)

        # remove false positive cross-junction concat_kmer
        # those not appear in the junction_kmer_df are all cross-junction kmer
        junction_kmer_df_all = get_unique_kmer(junction_file_name,zero_filter=False,cross_junction=False)
        final_junction_kmer_df = pd.concat((concat_kmer_df,junction_kmer_df_all))

        kmer_df_set = kmer_df_set.union(set(final_junction_kmer_df['kmer']))
    result_file = os.path.join(data_dir,'integrate_kmer_without_any_filter')
    f = open(result_file,'w')
    for kmer in kmer_df_set:
        f.write(kmer+'\n')
    return kmer_df_set

# create background kmer set
total_kmer_set = set()
TCGA_reference = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/ref'
Gtex_reference = '/cluster/work/grlab/projects/TCGA/immunopepper_gtex_rerun/TCGA-A2-A0D0'
ref_annotation = '/cluster/work/grlab/projects/TCGA/immunopepper_gtex_rerun/ref_annotation'

total_kmer_set = total_kmer_set.union(generate_unique_kmer_and_write(TCGA_reference,'ref'))
total_kmer_set = total_kmer_set.union(generate_unique_kmer_and_write(Gtex_reference,'ref'))
total_kmer_set = total_kmer_set.union(generate_unique_kmer_and_write(ref_annotation,'ref'))


result_dir = '/cluster/work/grlab/projects/TCGA/immunopepper_gtex_rerun'
all_dir = list(filter(lambda dir_name: dir_name.startswith('TCGA'),os.listdir(result_dir)))

for dir_name in all_dir:
    print(dir_name)
    data_dir = os.path.join(result_dir,dir_name)
    total_kmer_set = total_kmer_set.union(generate_unique_kmer_and_write(data_dir,'variant'))
    print_memory_diags()
    sys.stdout.flush()

result_file = os.path.join('/cluster/work/grlab/projects/TCGA/immunopepper_rerun/','integrate_kmer_without_any_filter')
f = open(result_file,'w')
for kmer in total_kmer_set:
    f.write(kmer+'\n')

