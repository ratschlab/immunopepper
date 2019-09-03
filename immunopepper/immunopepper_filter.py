"""
Apply different filter mechanism on kmer tsv files
"""
import pandas as pd
import os
def immunopepper_filter(arg):
    junction_kmer_tsv_path = arg.junction_kmer_tsv_path
    output_file_path = arg.output_file_path
    verbose = arg.verbose
    kmer_df = pd.read_csv(junction_kmer_tsv_path,sep='\t')
    if arg.cross_junction:
        if verbose > 1:
            print('apply cross junction filter')
        kmer_df = kmer_df[kmer_df['is_crossjunction']]

    if arg.seg_expr:
        seg_expr_thre = arg.seg_expr_thre
        if verbose > 1:
            print('apply segment expression filter, threshold is {}'.format(seg_expr_thre))
        kmer_df = kmer_df[kmer_df['seg_expr']>seg_expr_thre]

    if arg.junc_expr:
        # if we want to filter based on junction expression
        # we actually also do cross_junction filter because
        # only cross junction kmers have junction expression
        junc_expr_thre = arg.junc_expr_thre
        if verbose > 1:
            print('apply junction expression filter, threshold is {}'.format(junc_expr_thre))
        kmer_df = kmer_df[kmer_df['is_crossjunction']]
        kmer_df['junction_expr'] = pd.to_numeric(kmer_df['junction_expr'])
        kmer_df = kmer_df[kmer_df['junction_expr']>junc_expr_thre]
    kmer_df.to_csv(output_file_path,sep='\t',index=False)
    if verbose:
        print("Apply filter to {} and save result to {}".format(junction_kmer_tsv_path,output_file_path))




