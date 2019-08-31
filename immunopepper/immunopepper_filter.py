"""
Apply different filter mechanism on kmer tsv files
"""
import pandas as pd
import os
def immunopepper_filter(arg):
    kmer_tsv_path = arg.kmer_tsv_path
    output_dir = arg.output_dir
    verbose = arg.verbose
    item = kmer_tsv_path.split('/')
    kmer_csv_file_name = item[-1]
    if len(output_dir) == 0:
        output_dir = '/'.join(item[:-1])
    output_file_name = kmer_csv_file_name
    kmer_df = pd.read_csv(kmer_tsv_path,sep='\t')
    if arg.cross_junction:
        if verbose > 1:
            print('apply cross junction filter')
        kmer_df = kmer_df[kmer_df['is_crossjunction']]
        output_file_name = 'cj_'+output_file_name

    if arg.seg_expr:
        seg_expr_thre = arg.seg_expr_thre
        if verbose > 1:
            print('apply segment expression filter, threshold is {}'.format(seg_expr_thre))
        kmer_df = kmer_df[kmer_df['seg_expr']>seg_expr_thre]
        output_file_name = 'seg_expr_'+str(seg_expr_thre)+'_'+ output_file_name

    if arg.junc_expr:
        junc_expr_thre = arg.junc_expr_thre
        if verbose > 1:
            print('apply junction expression filter, threshold is {}'.format(junc_expr_thre))
        kmer_df['junction_expr'] = pd.to_numeric(kmer_df['junction_expr'])
        kmer_df = kmer_df[kmer_df['junction_expr']>junc_expr_thre]
        output_file_name = 'junc_expr_'+str(junc_expr_thre)+'_'+ output_file_name
    output_file_path = os.path.join(output_dir,output_file_name)
    kmer_df.to_csv(output_file_path,sep='\t',index=False)
    if verbose:
        print("Apply filter to {} and save result to {}".format(kmer_tsv_path,output_file_path))




