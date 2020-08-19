"""
Apply different filter mechanism on kmer tsv files
"""
import logging
import pandas as pd
import sys


from .constant import NOT_EXIST
from .io_ import save_pd_toparquet
from .io_ import read_pq_with_pd


def mode_filter(arg):
    logging.info(">>>>>>>>> filter: Start")
    junction_kmer_tsv_path = arg.junction_kmer_tsv_path
    output_file_path = arg.output_file_path
    if arg.infer_dna_pos:
        logging.info("infer exact dna pos for each kmer output so that we could apply rna-seq based filter")
        if not arg.meta_file_path:
            logging.error("no metadata file provided. Exit")
            sys.exit(-1)
        output_file_path = get_start_pos_for_kmer_unfiltered(arg.meta_file_path, junction_kmer_tsv_path, output_file_path,
                                          is_compressed=arg.compressed)
        logging.info("Infer dna pos for {} and save result to {}".format(junction_kmer_tsv_path,output_file_path))
        logging.info(">>>>>>>>> filter: Finish\n")
        return

    kmer_df = read_pq_with_pd(junction_kmer_tsv_path)



    if arg.cross_junction and kmer_df.shape[0]:
        logging.info('Apply cross junction filter')
        initial_shape = kmer_df.shape[0]
        is_cross_junction = kmer_df['is_cross_junction'].apply(lambda x: 'True' in x.split('/'))
        kmer_df = kmer_df[is_cross_junction]
        final_shape = kmer_df.shape[0]
        logging.info('...Initial size: {} kmers, Final size: {} kmers'.format(initial_shape, final_shape))

    if arg.seg_expr and kmer_df.shape[0]:
        seg_expr_thre = arg.seg_expr_thresh
        logging.info('Apply segment expression filter, threshold is {}'.format(seg_expr_thre))
        initial_shape = kmer_df.shape[0]
        seg_expr = kmer_df['segment_expr'].apply(lambda x: max(x.split('/')))
        seg_expr = seg_expr.astype('float')
        kmer_df = kmer_df[seg_expr>seg_expr_thre]
        final_shape = kmer_df.shape[0]
        logging.info('...Initial size: {} kmers, Final size: {} kmers'.format(initial_shape, final_shape))

    if arg.junc_expr and kmer_df.shape[0]:
        # if we want to filter based on junction expression
        # we actually also do cross_junction filter because
        # only cross junction kmers have junction expression
        junc_expr_thre = arg.junc_expr_thresh
        logging.info('Apply junction expression filter, threshold is {}'.format(junc_expr_thre))
        initial_shape = kmer_df.shape[0]
        is_cross_junction = kmer_df['is_cross_junction'].apply(lambda x: 'True' in x.split('/'))
        kmer_df = kmer_df[is_cross_junction]
        junction_expr = kmer_df['junction_expr'].apply(lambda x: max(x.split('/')))
        junction_expr = junction_expr.astype('float')
        kmer_df = kmer_df[junction_expr>junc_expr_thre]
        final_shape = kmer_df.shape[0]
        logging.info('...Initial size: {} kmers, Final size: {} kmers'.format(initial_shape, final_shape))

    if arg.meta_file_path:
        initial_shape = kmer_df.shape[0]
        meta_file_path = arg.meta_file_path
        meta_df = read_pq_with_pd(meta_file_path)
        meta_initial = meta_df.shape[0]
        total_keep_id = set(meta_df['id'])
        if arg.peptide_annotated:
            # annotated if 1 or 1/0 (metadata collapsed per gene)
            if arg.peptide_annotated in '1':
                keep_id = meta_df.loc[meta_df['peptide_annotated'].apply(lambda x: '1' in x.split('/')), 'id']

            else:
                keep_id = meta_df.loc[meta_df['peptide_annotated']=='0', 'id']
            total_keep_id = total_keep_id.intersection(keep_id)
            logging.info('Apply peptide_annotated filter, value is {}'.format(arg.peptide_annotated))
            logging.info('... {}/{} peptides satisfy condition'.format(len(keep_id), meta_initial))

        if arg.junction_annotated:
            # annotated if 1 or 1/0 (metadata collapsed per gene)
            if arg.junction_annotated in '1':
                keep_id = meta_df.loc[meta_df['junction_annotated'].apply(lambda x: '1' in x.split('/')), 'id']

            else:
                keep_id = meta_df.loc[meta_df['junction_annotated'] == '0', 'id']
            total_keep_id = total_keep_id.intersection(keep_id)
            logging.info('Apply junction_annotated filter, value is {}'.format(arg.junction_annotated))
            logging.info('... {}/{} peptides satisfy condition'.format(len(keep_id), meta_initial))

        if arg.has_stop_codon:
            # has_stop_codon if strictly 1 (metadata collapsed per gene)
            if arg.has_stop_codon in '1':
                keep_id = meta_df.loc[meta_df['has_stop_codon'] == '1', 'id']
            else:
                keep_id = meta_df.loc[meta_df['has_stop_codon'].apply(lambda x: '0' in x.split('/')), 'id']
            total_keep_id = total_keep_id.intersection(keep_id)
            logging.info('Apply has_stop_codon filter, value is {}'.format(arg.has_stop_codon))
            logging.info('... {}/{} peptides satisfy condition'.format(len(keep_id), meta_initial))

        if arg.is_in_junction_list:
            if '.' not in meta_df['is_in_junction_list'].values: # no whitelist provided
                # is_in_junction_list if 1 or 0/1 (metadata collapsed per gene)
                if arg.is_in_junction_list in '1':
                    keep_id = meta_df.loc[meta_df['is_in_junction_list'].apply(lambda x: '1' in x.split('/')), 'id']

                else:
                    keep_id = meta_df.loc[meta_df['is_in_junction_list'] == '0', 'id']
                total_keep_id = total_keep_id.intersection(keep_id)
                logging.info('Apply junction whitelist filter, value is {}'.format(arg.is_in_junction_list))
                logging.info('... {}/{} peptides satisfy condition'.format(len(keep_id), meta_initial))

        if arg.is_isolated:
            # is_isolated if strict 1 (metadata collapsed per gene)
            if arg.is_isolated in '1':
                keep_id = meta_df.loc[meta_df['is_isolated'] == '1', 'id']
            else:
                keep_id = meta_df.loc[meta_df['is_isolated'].apply(lambda x: '0' in x.split('/')), 'id']
            total_keep_id = total_keep_id.intersection(keep_id)
            logging.info('Apply is_isolated filter, value is {}'.format(arg.is_isolated))
            logging.info('... {}/{} peptides satisfy condition'.format(len(keep_id), meta_initial))

        #keep_id = {id for concat_id in total_keep_id for id in concat_id.split('/')}
        total_keep_id = '/'.join(total_keep_id)
        keep_kmer = {idx for idx, concat_id in enumerate(kmer_df['id']) for id in concat_id.split('/') if id in total_keep_id }
        logging.info('Final table {} kmers retained after filtering on metadata'.format(len(keep_kmer)))
        kmer_df = kmer_df.iloc[list(keep_kmer)]
        final_shape = kmer_df.shape[0]
        logging.info('...Initial size: {} kmers, Final size: {} kmers'.format(initial_shape, final_shape))

    if kmer_df.shape[0]:
        if arg.compressed:
            compression = 'gzip'
        else:
            compression = None

        save_pd_toparquet(output_file_path, kmer_df,
                      compression=compression, verbose=True)

    logging.info("Apply filter to {} and save result to {}".format(junction_kmer_tsv_path,output_file_path))
    logging.info(">>>>>>>>> filter: Finish\n")


def get_start_pos_for_kmer_unfiltered(meta_file_path,junction_kmer_tsv_path,output_file_path,is_compressed=True):
    """
    Given kmer_junction_file, return the exact dna positions that output the given kmer.
    Parameters
    ----------
    meta_file_path: _metadata.tsv(.gz) file outputed by build mode.
    junction_kmer_tsv_path: _junction_kmer.txt by build mode.
    output_file_path: the specified output file path.
    is_compressed: Bool. If output compressed file.

    Return:
    output_file_path: str. return the true output file name (in some cases will add .gz)
    """
    meta_df = pd.read_csv(meta_file_path,sep='\t',usecols=[0,3,4,11,13])
    if is_compressed:
        if not output_file_path.endswith('.gz'):
            output_file_path += '.gz'
    new_junction_file = gz_and_normal_open(output_file_path,'w')
    old_junction_file = gz_and_normal_open(junction_kmer_tsv_path,'r')
    headline = next(old_junction_file)
    new_headline = headline.strip()+'\texact_kmer_pos\n'
    new_junction_file.write(new_headline)
    count = 0
    prev_output_id = None
    line_buffer = []
    for line_id,line in enumerate(old_junction_file):
        if line_id % 100000 == 0:
            print(line_id)
        items = line.strip().split('\t')
        output_id = items[1]
        kmer_len = len(items[0])
        if prev_output_id != output_id:
            if line_buffer: # time to clear buffer
                deal_with_duplicate(line_buffer,num_dup,cur_modi_coord,vertex_len,strand,kmer_len,chr,cur_variant_comb,new_junction_file)
                line_buffer = []
            cur_meta_line = meta_df[meta_df['output_id'] == output_id]
            num_dup = len(cur_meta_line)
            if num_dup > 1: # initialize the buffer
                line_buffer.append(line)
            cur_modi_coord = cur_meta_line['modified_exons_coord'].values[0].split(';')
            cur_modi_coord = [int(coord) if coord != NOT_EXIST else coord for coord in cur_modi_coord]
            vertex_len = [int(cur_modi_coord[2*i+1])-int(cur_modi_coord[2*i]) if cur_modi_coord[2*i+1] != NOT_EXIST else 0 for i in range(len(cur_modi_coord)//2)]
            strand = cur_meta_line['gene_strand'].values[0]
            chr = cur_meta_line['gene_chr'].values[0]
            cur_variant_comb = cur_meta_line['variant_comb'].values[0]
            prev_output_id = output_id
            count = 0
        elif line_buffer: # the same id with previous, add to the buffer
            line_buffer.append(line)
            continue
        else: # no buffer initialized, use count directly
            count += 1
        pos_list = get_start_pos_from_count(count,cur_modi_coord,vertex_len,strand,kmer_len)
        pos_list = [str(pos) for pos in pos_list]
        exact_kmer_pos = str(chr)+'_'+strand+'_'+str(cur_variant_comb)+'_'+';'.join(pos_list)
        new_line = line.strip()+'\t'+exact_kmer_pos+'\n'
        new_junction_file.write(new_line)
    new_junction_file.close()
    return output_file_path


def deal_with_duplicate(line_list,k,cur_modi_coord,vertex_len,strand,kmer_len,chr,cur_variant_comb,new_junction_file):
    assert len(line_list) % k == 0
    uniq_line_len = len(line_list) // k
    for line_id,line in enumerate(line_list[:uniq_line_len]):
        pos_list = get_start_pos_from_count(line_id, cur_modi_coord, vertex_len, strand, kmer_len)
        pos_list = [str(pos) for pos in pos_list]
        exact_kmer_pos = str(chr) + '_' + strand + '_' + str(cur_variant_comb) + '_' + ';'.join(pos_list)
        new_line = line.strip() + '\t' + exact_kmer_pos + '\n'
        new_junction_file.write(new_line)


def get_start_pos_from_count(count, coord_list, vertex_len, strand, kmer_len):
    """
    From the kmer count if infer the start position and end position of that kmer
    Parameters
    ----------
    count: int. the *count*th kmer translated from the peptide.
    coord_list: list. coord list for the translated peptide.
    vertex_len: list. length of each vertex.
    strand: str. '+' or '-'.
    kmer_len: 3*k for k-mer

    Returns
    -------
    pos_list: list. start and end position of the dna that translates the given kmer. If
        spanning over k vertices, len(pos_list) = 2*k.

    """
    offset = count*3
    pos_list = []
    if strand == '+':
        cur_sum = vertex_len[0]
        i = 1
        while cur_sum < offset:
            cur_sum += vertex_len[i]
            i += 1
        start_pos = coord_list[2*(i-1)]+offset-sum(vertex_len[:i-1])
        pos_list.append(start_pos)
        end_pos = start_pos+kmer_len*3
        while end_pos > coord_list[2*(i-1)+1]:
            pos_list.append(coord_list[2*(i-1)+1]) # pos_list[-2], the stop pos of last vertex
            pos_list.append(coord_list[2*i]) # pos_list[-1], the start pos of cur vertex
            end_pos = pos_list[-1]+end_pos-pos_list[-2] #
            i += 1 # move the point to the next vertex
        pos_list.append(end_pos)
    else:
        cur_sum = vertex_len[0]
        i = 1
        while cur_sum < offset:
            cur_sum += vertex_len[i]
            i += 1
        start_pos = coord_list[2*(i-1)+1]-(offset-sum(vertex_len[:i-1]))
        pos_list.append(start_pos)
        end_pos = start_pos-kmer_len*3
        while end_pos < coord_list[2*(i-1)]:
            pos_list.append(coord_list[2*(i-1)]) # pos_list[-2], the stop pos of last vertex
            pos_list.append(coord_list[2*i+1]) # pos_list[-1], the start pos of cur vertex
            end_pos = pos_list[-1]-(pos_list[-2]-end_pos) #
            i += 1 # move the point to the next vertex
        pos_list.append(end_pos)
    return pos_list
