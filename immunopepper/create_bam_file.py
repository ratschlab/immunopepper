from immuno_preprocess import preprocess_ann
from utils import complementary_seq
import numpy as np
import os


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

    new_dna_list = dna_list+list(complementary_seq(dna_list[::-1])) # create the same
    dna = ''.join(new_dna_list)
    seq_file = os.path.join(data_dir,'test1.fa')
    f = open(seq_file,'w')
    f.write('>X dna:chromosome chromosome:GRCh37:X:1:155270560:1'+'\n')
    f.write(dna)
    f.close()
    return dna


def create_full_from_pos(data_dir, L):
    posgtf_file = os.path.join(data_dir, 'test1pos.gtf')
    posmaf_file = os.path.join(data_dir, 'test1pos.maf')
    gtf_file = os.path.join(data_dir, 'test1.gtf')
    maf_file = os.path.join(data_dir, 'test1.maf')
    f = open(posgtf_file,'r')
    lines = f.readlines()
    new_line_list = []
    for line in lines:
        item = line.split('\t')
        start_pos = 2*L+1-int(item[3])
        end_pos = 2*L+1-int(item[4])
        new_detail = item[8].replace('GENE1', 'GENE2')
        new_detail = new_detail.replace('TRANS1', 'TRANS2')
        item[3] = str(end_pos)
        item[4] = str(start_pos)
        item[8] = new_detail
        item[6] = '-'
        new_line = '\t'.join(item)
        new_line_list.append(new_line)
    f = open(gtf_file,'w')
    f.writelines(lines+new_line_list)
    f.close()

    f = open(posmaf_file, 'r')
    lines = f.readlines()
    new_line_list = []
    for line in lines[1:]:
        item = line.split('\t')
        pos = item[5]
        pre_base = item[10]
        new_base = item[12]
        neg_pre_base = complementary_seq(pre_base)
        neg_new_base = complementary_seq(new_base)
        new_pos = 2*L+1-int(pos)
        item[5] = str(new_pos)
        item[6] = str(new_pos)
        item[7] = '-'
        item[10] = neg_pre_base
        item[11] = neg_pre_base
        item[12] = neg_new_base
        new_line = '\t'.join(item)
        new_line_list.append(new_line)
    f = open(maf_file,'w')
    f.writelines(lines+new_line_list)
    f.close()
    #f = open(posvcf_file, 'r')
    #todo complete the pos maf part

def create_bam_file(data_dir,seed=1):
    np.random.seed(seed)
    ref_path = os.path.join(data_dir,'test1.fa')
    f = open(ref_path,'r')
    ref_seq = f.readlines()[1]
    ann_path = os.path.join(data_dir, 'test1.gtf')
    gene_cds_begin_dict, gene_to_transcript_table, transcript_to_cds_table = preprocess_ann(ann_path)

    expr = [10, 0, 25, 15, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 20]
    reads_L = 15
    cover = 100

    total_reads_list = []
    for i,key_value_pair in enumerate(transcript_to_cds_table.items()):
        if i < 6:
            exon_coord = key_value_pair[1]
            ts = ''
            for coord_pair in exon_coord:
                ts += ref_seq[coord_pair[0]+coord_pair[2]:coord_pair[1]+1]
            N = expr[i]*cover
            pos = np.random.randint(0,len(ts)-reads_L,N)
            reads_list = [ts[ipos:ipos+reads_L] for ipos in pos]
            total_reads_list.extend(reads_list)
    output_dir = os.path.join(data_dir,'test1_1.fq')
    f = open(output_dir,'w')
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

data_dir = '/Users/jiayu/PycharmProjects/CBM_RA/projects2018_immunopepper/tests/data'
create_test_genome(150, [73, 102],data_dir)
print("Successfully create test genome")
create_full_from_pos(data_dir, 150)
print("Successfully create gtf, maf and vcf file")
create_bam_file(data_dir)
print("Successfully create bam segments")
