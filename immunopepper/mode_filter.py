"""
Apply different filter mechanism on kmer tsv files
"""
import logging
import sys

import pandas as pd

from .constant import NOT_EXIST
from .io_ import gz_and_normal_open

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

    kmer_df = pd.read_csv(junction_kmer_tsv_path,sep='\t')
    if arg.cross_junction:
        logging.info('apply cross junction filter')
        kmer_df = kmer_df[kmer_df['is_crossjunction']]

    if arg.seg_expr:
        seg_expr_thre = arg.seg_expr_thresh
        logging.info('apply segment expression filter, threshold is {}'.format(seg_expr_thre))
        kmer_df = kmer_df[kmer_df['seg_expr']>seg_expr_thre]

    if arg.junc_expr:
        # if we want to filter based on junction expression
        # we actually also do cross_junction filter because
        # only cross junction kmers have junction expression
        junc_expr_thre = arg.junc_expr_thresh
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
            logging.info('apply peptide_annotated filter, value is {}'.format(arg.peptide_annotated))

        if arg.junction_annotated:
            if int(arg.junction_annotated):
                keep_id = meta_df[meta_df['junction_annotated'].isin(['1'])]['output_id']
            else:
                keep_id = meta_df[meta_df['junction_annotated'].isin(['0'])]['output_id']
            total_keep_id = total_keep_id.intersection(keep_id)
            logging.info('apply junction_annotated filter, value is {}'.format(arg.junction_annotated))

        if arg.has_stop_codon:
            keep_id = meta_df[meta_df['has_stop_codon']==int(arg.has_stop_codon)]['output_id']
            total_keep_id = total_keep_id.intersection(keep_id)
            logging.info('apply has_stop_codon filter, value is {}'.format(arg.has_stop_codon))

        if arg.is_in_junction_list:
            if int(arg.junction_annotated):
                keep_id = meta_df[meta_df['junction_annotated'].isin(['1'])]['output_id']
            else:
                keep_id = meta_df[meta_df['junction_annotated'].isin(['0'])]['output_id']
            total_keep_id = total_keep_id.intersection(keep_id)
            logging.info('apply junction whitelist filter, value is {}'.format(arg.is_in_junction_list))

        if arg.is_isolated:
            keep_id = meta_df[meta_df['is_isolated']==int(arg.is_isolated)]['output_id']
            total_keep_id = total_keep_id.intersection(keep_id)
            logging.info('apply is_isolated filter, value is {}'.format(arg.is_isolated))

        kmer_df = kmer_df[kmer_df['gene_name'].isin(total_keep_id)]
    if arg.compressed:
        kmer_df.to_csv(output_file_path, sep='\t', index=False,compression='gzip')
    else:
        kmer_df.to_csv(output_file_path,sep='\t',index=False)

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
