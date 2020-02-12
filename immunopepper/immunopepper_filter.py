"""
Apply different filter mechanism on kmer tsv files
"""
import pandas as pd
import logging

def immunopepper_filter(arg):
    logging.info(">>>>>>>>> filter: Start")
    junction_kmer_tsv_path = arg.junction_kmer_tsv_path
    output_file_path = arg.output_file_path
    verbose = arg.verbose
    kmer_df = pd.read_csv(junction_kmer_tsv_path,sep='\t')
    if arg.cross_junction:
        if verbose > 1:
            logging.info('apply cross junction filter')
        kmer_df = kmer_df[kmer_df['is_crossjunction']]

    if arg.seg_expr:
        seg_expr_thre = arg.seg_expr_thresh
        if verbose > 1:
            logging.info('apply segment expression filter, threshold is {}'.format(seg_expr_thre))
        kmer_df = kmer_df[kmer_df['seg_expr']>seg_expr_thre]

    if arg.junc_expr:
        # if we want to filter based on junction expression
        # we actually also do cross_junction filter because
        # only cross junction kmers have junction expression
        junc_expr_thre = arg.junc_expr_thresh
        if verbose > 1:
            logging.info('apply junction expression filter, threshold is {}'.format(junc_expr_thre))
        kmer_df = kmer_df[kmer_df['is_crossjunction']]
        kmer_df['junction_expr'] = pd.to_numeric(kmer_df['junction_expr'])
        kmer_df = kmer_df[kmer_df['junction_expr']>junc_expr_thre]

    if arg.meta_file_path:
        meta_file_path = arg.meta_file_path
        meta_df = pd.read_csv(meta_file_path,sep='\t')
        total_keep_id = set(meta_df['output_id'])

        if arg.peptide_annotated:
            keep_id = meta_df[meta_df['peptide_annotated']==int(arg.peptide_annotated)]['output_id']
            total_keep_id = total_keep_id.intersection(keep_id)
            if verbose > 1:
                logging.info('apply peptide_annotated filter, value is {}'.format(arg.peptide_annotated))

        if arg.junction_annotated:
            if int(arg.junction_annotated):
                keep_id = meta_df[meta_df['junction_annotated'].isin(['1','1;0','0;1','1;1'])]['output_id']
            else:
                keep_id = meta_df[meta_df['junction_annotated'].isin(['0','0;0'])]['output_id']
            total_keep_id = total_keep_id.intersection(keep_id)
            if verbose > 1:
                logging.info('apply junction_annotated filter, value is {}'.format(arg.junction_annotated))

        if arg.has_stop_codon:
            keep_id = meta_df[meta_df['has_stop_codon']==int(arg.has_stop_codon)]['output_id']
            total_keep_id = total_keep_id.intersection(keep_id)
            if verbose > 1:
                logging.info('apply has_stop_codon filter, value is {}'.format(arg.has_stop_codon))

        if arg.is_in_junction_list:
            keep_id = meta_df[meta_df['is_in_junction_list']==int(arg.is_in_junction_list)]['output_id']
            total_keep_id = total_keep_id.intersection(keep_id)
            if verbose > 1:
                logging.info('apply junction whitelist filter, value is {}'.format(arg.is_in_junction_list))

        if arg.is_isolated:
            keep_id = meta_df[meta_df['is_isolated']==int(arg.is_isolated)]['output_id']
            total_keep_id = total_keep_id.intersection(keep_id)
            if verbose > 1:
                logging.info('apply is_isolated filter, value is {}'.format(arg.is_isolated))

        kmer_df = kmer_df[kmer_df['gene_name'].isin(total_keep_id)]
    if arg.compressed:
        kmer_df.to_csv(output_file_path, sep='\t', index=False,compression='gzip')
    else:
        kmer_df.to_csv(output_file_path,sep='\t',index=False)

    if verbose:
        logging.info("Apply filter to {} and save result to {}".format(junction_kmer_tsv_path,output_file_path))
    logging.info(">>>>>>>>> filter: Finish\n")



