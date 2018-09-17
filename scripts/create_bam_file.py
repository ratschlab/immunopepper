from immunopepper.immuno_preprocess import preprocess_ann
from immunopepper.utils import complementary_seq
import numpy as np
import os
from collections import Counter

# steps to add new transcript
# 1. modify test1pos.gtf file
# 2. modify create_test_genome file to make tailor-made gene
# 3. delete quick_test_data/spladder
# 4. create new splicegraph with command tools
# 4. run main_immuno.py
def create_test_genome(L,stop_position,data_dir,seed=1):
    np.random.seed(seed)
    number_list = np.random.randint(0,4,L)
    map_dict = {0:'A', 1:'G', 2:'C', 3:'T'}
    dna_list = [map_dict[i] for i in number_list]
    for pos in stop_position:
        dna_list[pos:pos+3] = 'TAG'

    # modify
    dna_list[69] = 'C'  # remove stop codon
    dna_list[98] = 'G'  # remove stop codon
    dna_list[39] = 'T'  # create stop codon in mutation mode

    # create GT/AG for splicing site
    dna_list[25] = 'G'  # in case the dna_list[24:27] = 'TAG'
    dna_list[26:28] = 'GT'
    dna_list[29:31] = 'GT'
    dna_list[36:38] = 'AG'
    dna_list[50:52] = 'GT'
    dna_list[59:61] = 'AG'
    dna_list[64:66] = 'AG'
    dna_list[75:77] = 'GT'
    dna_list[85:87] = 'AG'

    pos_dna = ''.join(dna_list)
    neg_dna = complementary_seq(pos_dna[::-1])
    pos_seq_file = os.path.join(data_dir,'test1pos.fa')
    neg_seq_file = os.path.join(data_dir,'test1neg.fa')
    f_pos = open(pos_seq_file, 'w')
    f_neg = open(neg_seq_file, 'w')
    f_pos.write('>X dna:chromosome chromosome:GRCh37:X:1:155270560:1'+'\n')
    f_pos.write(pos_dna)
    f_pos.close()
    f_neg.write('>X dna:chromosome chromosome:GRCh37:1:1:155270560:1'+'\n')
    f_neg.write(neg_dna)
    f_neg.close()
    return pos_dna


def create_full_from_pos(data_dir, L):
    posgtf_file = os.path.join(data_dir, 'test1pos.gtf')
    posmaf_file = os.path.join(data_dir, 'test1pos.maf')
    posvcf_file = os.path.join(data_dir, 'test1pos.vcf')
    neggtf_file = os.path.join(data_dir, 'test1neg.gtf')
    negmaf_file = os.path.join(data_dir, 'test1neg.maf')
    negvcf_file = os.path.join(data_dir, 'test1neg.vcf')

    f = open(posgtf_file, 'r')
    lines = f.readlines()
    new_line_list = []
    for line in lines:
        item = line.split('\t')
        start_pos = L+1-int(item[3])
        end_pos = L+1-int(item[4])
        #new_detail = item[8].replace('GENE1', 'GENE2')
        #new_detail = new_detail.replace('TRANS1', 'TRANS2')
        item[3] = str(end_pos)
        item[4] = str(start_pos)
        item[6] = '-'
        new_line = '\t'.join(item)
        new_line_list.append(new_line)
    f = open(neggtf_file,'w')
    f.writelines(new_line_list)
    f.close()

    f = open(posmaf_file, 'r')
    lines = f.readlines()
    new_line_list = [lines[0]]
    for line in lines[1:]:
        item = line.split('\t')
        pos = item[5]
        ref_base = item[10]
        mut_base = item[12]
        sample = item[15]
        neg_pre_base = complementary_seq(ref_base)
        neg_new_base = complementary_seq(mut_base)
        new_sample = sample.replace('pos', 'neg')
        new_pos = L+1-int(pos)
        item[5] = str(new_pos)
        item[6] = str(new_pos)
        item[7] = '-'
        item[10] = neg_pre_base
        item[11], item[12] = neg_pre_base, neg_new_base
        item[15], item[16] = new_sample, new_sample
        new_line = '\t'.join(item)
        new_line_list.append(new_line)
    f = open(negmaf_file,'w')
    f.writelines(new_line_list)
    f.close()

    f = open(posvcf_file, 'r')
    lines = f.readlines()
    new_line_list = lines[:29]
    headline = lines[29].replace('pos','neg')
    new_line_list.extend(headline)
    for line in lines[30:]:
        item = line.split('\t')
        pos = item[1]
        ref_base = item[3]
        mut_base = item[4]
        neg_pre_base = complementary_seq(ref_base)
        neg_new_base = complementary_seq(mut_base)
        neg_pos = L+1-int(pos)
        item[1] = str(neg_pos)
        item[3] = neg_pre_base
        item[4] = neg_new_base
        new_line = '\t'.join(item)
        new_line_list.extend(new_line)
    f = open(negvcf_file,'w')
    f.writelines(new_line_list)
    f.close()


def create_bam_file(data_dir, seed=1):
    np.random.seed(seed)
    posref_path = os.path.join(data_dir, 'test1pos.fa')
    with open(posref_path, 'r') as fpos:
        posref_seq = fpos.readlines()[1]
    posann_path = os.path.join(data_dir, 'test1pos.gtf')
    pos_genetable = preprocess_ann(posann_path)

    expr = [10, 0, 25, 15, 10]
    reads_L = 15
    cover = 100

    expr_file = os.path.join(data_dir, 'expr_ts.txt')
    f_expr = open(expr_file, 'w')
    total_reads_list = []
    ts_dict = {}
    expr_list = []
    for i, key_value_pair in enumerate(pos_genetable.ts_to_cds.items()):
        ts_name = key_value_pair[0]
        pos_exon_coord = key_value_pair[1]
        ts = ''
        for coord_pair in pos_exon_coord:
            start_pos = coord_pair[0]
            read_frame = coord_pair[2]
            stop_pos = coord_pair[1]
            ts += posref_seq[start_pos+read_frame:stop_pos+1]
        N = expr[i]*cover
        pos = np.random.randint(0,len(ts)-reads_L,N)
        pos_counter = Counter(pos)
        expr_arr = np.zeros(len(ts))
        expr_list.append(expr_arr)
        ts_dict[ts_name] = ts
        for ipos in pos_counter.keys():
            expr_arr[ipos:ipos+reads_L] += pos_counter[ipos]
        expr_arr = [str(int(i_expr)) for i_expr in expr_arr]
        new_line = ts_name+'\t'+'\t'.join(expr_arr)+'\n'
        f_expr.write(new_line)
        reads_list = [ts[ipos:ipos+reads_L] for ipos in pos]
        total_reads_list.extend(reads_list)
    f_expr.close()
    output_dir = os.path.join(data_dir, 'test1_1.fq')
    f = open(output_dir, 'w')
    total_num = len(total_reads_list)
    for i in range(total_num-1, -1, -1):
        line1 = '@test1_' + str((i+1)*2) + '/1'+'\n'
        f.write(line1)
        line2 = total_reads_list[i]+'\n'
        f.write(line2)
        line3 = '+\n'
        f.write(line3)
        line4 = reads_L*'G'+'\n'
        f.write(line4)


data_dir = '/Users/jiayu/PycharmProjects/CBM_RA/projects2018_immunopepper/tests/test1/data'
create_test_genome(150, [73, 102], data_dir)
print("Successfully create test genome")
create_full_from_pos(data_dir, 150)
print("Successfully create gtf, maf and vcf file")
create_bam_file(data_dir)
print("Successfully create bam segments")
