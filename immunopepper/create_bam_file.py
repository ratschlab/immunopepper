from immuno_preprocess import preprocess_ann
from immuno_filter import find_background_transcript
import cPickle
import numpy as np

def create_bam_file():
    # get 4 transcripts
    # setting the coverage and read length
    # setting the reletive expression data
    # get the N
    ref_seq = 'GTAATGTGTAAGATGACGCACGCATACAGCCATTGGTCATGGGTAGCGGAAGAAGTTCGACTTCGGGTACGTATAGGTGCTAAGTCTTCGACGTTAGGGTGGTAGACCTTTTGTTACCCAATTGGGTAAGGCAGCCATGTTGATATACATATGTATATCAACATGGCTGCCTTACCCAATTGGGTAACAAAAGGTCTACCACCCTAACGTCGAAGACTTAGCACCTATACGTACCCGAAGTCGAACTTCTTCCGCTACCCATGACCAATGGCTGTATGCGTGCGTCATCTTACACATTAC'
    ann_path = 'quick_test_data/test1.gtf'
    gene_cds_begin_dict, gene_to_transcript_table, transcript_to_cds_table = preprocess_ann(ann_path)
    splice_path = 'quick_test_data/spladder/genes_graph_conf3.merge_graphs.pickle'
    (graph_data, graph_meta) = cPickle.load(open(splice_path,'r'))
    ts_list = find_background_transcript(graph_data[0],ref_seq,gene_to_transcript_table,transcript_to_cds_table)
    expr = [10, 0, 25, 15]
    reads_L = 15
    cover = 100
    total_reads_list = []
    for i,ts in enumerate(ts_list):
        N = expr[i]*cover
        pos = np.random.randint(0,len(ts)-reads_L,N)
        reads_list = [ts[ipos:ipos+reads_L] for ipos in pos]
        total_reads_list.extend(reads_list)
    f = open('test1_5.fq','w')
    total_num = len(total_reads_list)
    for i in range(total_num-1,-1,-1):
        line1 = '@test1_' + str((i+1)*2) + '/1'+'\n'
        f.write(line1)
        line2 = total_reads_list[i]+'\n'
        f.write(line2)
        line3 = '+\n'
        f.write(line3)
        line4 = reads_L*'G'+'\n'
        f.write(line4)
create_bam_file()