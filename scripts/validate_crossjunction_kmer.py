# validation work in March 9th
# validate 03/07/2019
import numpy as np
import gzip
from collections import Counter
sample_name = 'TCGA-13-1489'
mutation_mode = 'germline'

immunopepper_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/{}/{}_junction_kmer.txt'.format(sample_name,mutation_mode)
nora_result_file = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/analysis_pancan/' \
                   'ccell_rerun_2018/output/peptides.clean/split/cj_kmers/{}.cj.{}.cj_kmers_9.fa'.format(sample_name,mutation_mode)
meta_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/{}/{}_metadata.tsv'.format(sample_name,mutation_mode)
meta_gz_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/{}/{}.metadata.tsv.gz'.format(sample_name,mutation_mode)

def get_kmer_list(_str,k):
    if len(_str) < k:
        return [_str]
    else:
        kmer_list = [_str[i:i+k] for i in range(0,max(0,len(_str)-k+1))]
        return kmer_list


# parse nora's result
nora_pep_lines = open(nora_result_file,'r').readlines()
i = 0
pep_dict = {}
nora_kmer_list = []
while i < len(nora_pep_lines):
    if i % 100000 == 0:
        print(i)
    line = nora_pep_lines[i]
    gene_name = line.strip().split('_')[2]
    i += 1
    pep = nora_pep_lines[i].strip()
    nora_kmer_list.extend(get_kmer_list(pep,9))
    i += 1
nora_kmer_set = set(nora_kmer_list)

# construct immuno_pepper dict, only focus on the cross junction part
# kmer: -> immuno_pepper name
immunopepper_dict = {line.split('\t')[0]:line.split('\t')[1] for line in open(immunopepper_file,'r') if line.strip().split('\t')[-1] == 'True'}

# get all kmer returned by immunopepper
immunopepper_kmer = set(immunopepper_dict.keys())
common_kmer_set = nora_kmer_set.intersection(immunopepper_kmer)
additional_kmer_list = list(immunopepper_kmer.difference(nora_kmer_set))
miss_kmer_list = list(nora_kmer_set.difference(immunopepper_kmer))
s_summary = "{} common kmer, {} additional kmer and {} miss kmer".format(len(common_kmer_set),len(additional_kmer_list),len(miss_kmer_list))
print(s_summary)
# 59680/60457; 59680/111186 (germline)
# (somatic)

# check missing kmer
# germline mutation occur at the extrapolation part
# so the whole peptide is kept alghough they are the same with ref peptide near the cross junction place
immunopepper_ref_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/TCGA-13-1489/ref_junction_kmer.txt'
ref_kmer_dict = {line.split('\t')[0]:line.split('\t')[1] for line in open(immunopepper_ref_file,'r') if line.strip().split('\t')[-1] == 'True'}
ref_kmer_set = set(ref_kmer_dict.keys())
s_explain_missing = "{} kmers can not be found in reference kmer".format(len(set(miss_kmer_list).intersection(ref_kmer_set)))
print(s_explain_missing)
problem_id_set = list(set([immunopepper_dict[_kmer] for _kmer in additional_kmer_list]))

"""
will remove this part
"""
# explain the additioanal kmers
# get some index information (gene_name+vertex_id) from immunopepper dictionary
f_meta = open(meta_file,'r')
meta_dict = {}
for i,line in enumerate(f_meta):
    if i ==0:
        continue
    items = line.strip().split('\t')
    gene_name = items[2]
    vertex_id_list = items[15].split(',')
    key = gene_name+'_'+vertex_id_list[0]+'_'+vertex_id_list[1]
    if key in problem_id_set:
        stop_flag = int(items[9])
        isolated_flag = int(items[11])
        coord_list = items[14].split(';')
        start_v1 = coord_list[0]
        stop_v1 = coord_list[1]
        is_less_than3_flag = int((int(stop_v1)-int(start_v1)) < 3)
        if key not in meta_dict:
            meta_dict[key] = [(stop_flag,isolated_flag,is_less_than3_flag)]
        else:
            meta_dict[key].append((stop_flag,isolated_flag,is_less_than3_flag))

flag_tuples = list(meta_dict.values())
flag_array = np.array(flag_tuples,dtype=int)
flag_result = np.sum(flag_array,axis=0)
flag_individual = np.sum(flag_array,axis=1)
s_explain_additional = "the number of additional cases that can be explained by stop codon is {}, by isolated is {}, by short_than_3 is {}".format(flag_result[0],flag_result[1],flag_result[2])
print(s_explain_additional)

bad_id_array = np.where(flag_individual < 1)[0]
bad_gene_name_list = np.array(meta_dict.keys())[bad_id_array]
# cluster bad_result with gene
bad_gene_dict = {}
for gene_name_str in bad_gene_name_list:
    gene_name = gene_name_str.split('_')[0]
    if gene_name not in bad_gene_dict:
        bad_gene_dict[gene_name] = [gene_name_str]
    else:
        bad_gene_dict[gene_name].append(gene_name_str)
"""
"""

# load graph
import cPickle
anno_pickle_path = open('/cluster/work/grlab/projects/immunopepper/annotation_preprop.pickle','r')
graph,cds_dict = cPickle.load(anno_pickle_path)

# load genetable
from immunopepper.immuno_preprocess import preprocess_ann
ann_path = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/annotation/gencode.v19.annotation.hs37d5_chr.gff'
genetable = preprocess_ann(ann_path)

gene_id_dict = {}
for i,gene in enumerate(graph):
    gene_id_dict[gene.name] = i

# load sequence file
import Bio.SeqIO as BioIO
ref_path = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/sequence/genome.fa'
seq_dict = {}
interesting_chr = map(str, range(1, 23)) + ["X", "Y", "MT"]
print('Parsing genome sequence ...')
for record in BioIO.parse(ref_path, "fasta"):
    if record.id in interesting_chr:
        seq_dict[record.id] = str(record.seq).strip()

from immunopepper.utils import translate_dna_to_peptide
from immunopepper.utils import complementary_seq


# example 1
seq7 = seq_dict['7']
print(translate_dna_to_peptide(complementary_seq(seq7[2854510:2854891][::-1]))) # ('LSSLLDCRSADWLPP', True) Matthias implementation (start from extrapolation)
print(translate_dna_to_peptide(complementary_seq(seq7[2854510:2854564][::-1]))) # ('MKRRMFPRPCLARMPGSR', False) immunopepper implementation (start from cds)

# example 2
seq = seq_dict['6']
print(translate_dna_to_peptide(complementary_seq(seq[79787745:79787953][::-1]))) # ('LSSLLDCRSADWLPP', True) Matthias implementation (start from extrapolation)
print(translate_dna_to_peptide(complementary_seq(seq[79787745:79787785][::-1]))) # ('MKRRMFPRPCLARMPGSR', False) immunopepper implementation (start from cds)

seq = seq_dict['10']
print(translate_dna_to_peptide(seq[133918174:133918447])) # ('MKRRMFPRPCLARMPGSR', False) immunopepper implementation (start from cds)
print(translate_dna_to_peptide(seq[133918312:133918447])) # ('LSSLLDCRSADWLPP', True) Matthias implementation (start from extrapolation)

seq = seq_dict['17']
print(translate_dna_to_peptide(seq[48539544:48539663])) # ('MKRRMFPRPCLARMPGSR', False) immunopepper implementation (start from cds)
print(translate_dna_to_peptide(seq[48539546:48539663])) # ('LSSLLDCRSADWLPP', True) Matthias implementation (start from extrapolation)

# example 3 (one of the special 58)
seq = seq_dict['8']
print(translate_dna_to_peptide(complementary_seq(seq[101724589:101724685][::-1]+seq[101721686:101721947][::-1]))) # ('LSSLLDCRSADWLPP', True) Matthias implementation (start from extrapolation)
print(translate_dna_to_peptide(complementary_seq(seq[101724589:101724685][::-1]))) # ('LSSLLDCRSADWLPP', True) Matthias implementation (start from extrapolation)

# try to explain the left
explain_flag_list = []
unexplained_result = bad_gene_dict.values()
for gene_name_str_list in unexplained_result:
    if len(gene_name_str_list) == 1:
        items = gene_name_str_list[0].split('_')
        gene_name = items[0]
        first_v_id = int(items[1])
        second_v_id = int(items[2])
        gene = graph[gene_id_dict[gene_name]]
        strand = gene.strand
        v_coord = gene.splicegraph.vertices[:,first_v_id]
        orig_len = v_coord[1]-v_coord[0]
        additional_base = orig_len % 3
        seq = seq_dict[gene.chr]
        if strand == '+':
            has_stop_flag = translate_dna_to_peptide(seq[v_coord[0]+additional_base:v_coord[1]])[1]
        else:
            has_stop_flag = translate_dna_to_peptide(complementary_seq(seq[v_coord[0]:v_coord[1]-additional_base][::-1]))[1]
        explain_flag_list.append(int(has_stop_flag))
    else:
        explain_flag_list.append(-1)
print("{} cases can be explained by the stop codon introduced by extrapolation".format(np.sum(np.array(explain_flag_list)>0)))


# final validation:
# build the missing table via reference comparing
f_gt = open('/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/analysis_pancan/ccell_rerun_2018/output/peptides.clean/split/REFERENCE.cj.fa','r')
def get_matthias_ref_junction_dict(f_gt):
    # reconstrcut gt_dict. using four coords as the key
    # calculate the distribution of ball-bin
    gt_lines = f_gt.readlines()
    ref_junction_dict = {}
    mat_gene_coord_dict = {}
    i = 0
    while i < len(gt_lines)-1:
        headline = gt_lines[i].strip()
        items = headline.split('_')
        gene_name = items[2]
        coord_str_tuple = tuple([int(coord)-1 if a%2 ==0 else int(coord) for a, coord in enumerate(items[-4:])]) # remove coord correction
        i += 1
        peptide = gt_lines[i].strip()
        ref_junction_dict[(gene_name,coord_str_tuple)] = peptide
        if gene_name not in mat_gene_coord_dict:
            mat_gene_coord_dict[gene_name] = [coord_str_tuple]
        else:
            mat_gene_coord_dict[gene_name].append(coord_str_tuple)
        i += 1
    return ref_junction_dict,mat_gene_coord_dict

meta_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/TCGA-13-1489/ref_metadata.tsv'
def get_immunopepper_meta_dict(meta_file):
    f_meta = open(meta_file, 'r')
    f_meta.readline()
    meta_flag_dict = {}
    imm_gene_coord_dict = {}
    for line in f_meta:
        items = line.strip().split('\t')
        gene_name = items[2]
        stop_flag = int(items[9])
        isolated_flag = int(items[11])
        coord_str_tuple = tuple([int(coord) if coord != NOT_EXIST else NOT_EXIST for coord in items[14].split(';') ])
        vertex_id = items[15]
        key = gene_name + '_' + '_'.join(items[14].split(';'))+'_'+vertex_id
        start_v1 = coord_str_tuple[0]
        stop_v1 = coord_str_tuple[1]
        is_less_than3_flag = int((int(stop_v1) - int(start_v1)) < 3)
        meta_flag_dict[key] = (stop_flag, isolated_flag, is_less_than3_flag)
        if gene_name not in imm_gene_coord_dict:
            imm_gene_coord_dict[gene_name] = {(coord_str_tuple,vertex_id)}
        else:
            imm_gene_coord_dict[gene_name].add((coord_str_tuple,vertex_id))
    return meta_flag_dict,imm_gene_coord_dict


from immunopepper.constant import NOT_EXIST
def find_diff_coord(imm_gene_coord_dict,mat_gene_coord_dict):
    def fuzzy_comp_str2str(query_coord,another_coord):
        return np.sum(np.array(query_coord)-np.array(another_coord)) == 0

    def fuzzy_comp_str2list(query_coord,another_coord_list):
        if query_coord[2] == NOT_EXIST or len(another_coord_list) == 0:
            return ('',False)

        possible_match_id_list = np.where(np.array(another_coord_list)[:,0]==query_coord[0])[0]
        for possible_match_id in possible_match_id_list:
            coord_str_tuple = another_coord_list[possible_match_id]
            if coord_str_tuple[3] != NOT_EXIST:
                match = fuzzy_comp_str2str(query_coord,coord_str_tuple)
                if match:
                    return (coord_str_tuple,True)
        return ('',False)

    unique_imm_coord_dict = {}
    for i,gene_name in enumerate(imm_gene_coord_dict):
        if i % 1000 == 0:
            print(i)
        imm_coord_set = imm_gene_coord_dict[gene_name]
        if gene_name not in mat_gene_coord_dict:
            mat_coord_set = {}
        else:
            mat_coord_set = mat_gene_coord_dict[gene_name]
        for coord,vertex_id in imm_coord_set:
            similar_result,found = fuzzy_comp_str2list(coord,mat_coord_set)
            if not found:
                gene_name_coord_str = gene_name+'_'+'_'.join([str(icoord) for icoord in coord])+'_'+vertex_id
                if gene_name not in unique_imm_coord_dict:
                    unique_imm_coord_dict[gene_name] = [gene_name_coord_str]
                else:
                    unique_imm_coord_dict[gene_name].append(gene_name_coord_str)
    return unique_imm_coord_dict

meta_flag_dict, imm_gene_coord_dict = get_immunopepper_meta_dict(meta_file)
ref_junction_dict, mat_gene_coord_dict = get_matthias_ref_junction_dict(f_gt)
unique_imm_coord_dict = find_diff_coord(imm_gene_coord_dict,mat_gene_coord_dict)

# do a simple analysis of the result
miss_junc_pair_num = np.sum([len(item) for item in unique_imm_coord_dict.values()]) # 1103270
total_junc_pair_num = np.sum([len(item) for item in imm_gene_coord_dict.values()]) # 1622126
total_mat_junc_pair_num = np.sum([len(item) for item in mat_gene_coord_dict.values()]) # 5255539
print("{} missing junc pair in {} immuno total junc pairs and {} mat junc pairs".format(miss_junc_pair_num,total_junc_pair_num,total_mat_junc_pair_num))

# see detail information
unique_imm_flag_dict = {}
unique_flag_list = []
unique_gene_name_list = []
for gene_name, miss_junc_list in unique_imm_coord_dict.items():
    unique_imm_flag_dict[gene_name] = [(miss_junc,meta_flag_dict[miss_junc]) for miss_junc in miss_junc_list]
    unique_flag_list.extend([meta_flag_dict[miss_junc] for miss_junc in miss_junc_list])
    for miss_junc in miss_junc_list:
        items = miss_junc.split('_')
        vertex_id = items[-1].split(',')
        gene_name = items[0]
        gene_name_str = '_'.join((gene_name,vertex_id[0],vertex_id[1]))
        unique_gene_name_list.append(gene_name_str)

unique_flag_array = np.array(unique_flag_list)
flag_explain_result = np.sum(unique_flag_array, axis=0)
print("stop codon constributes to {}, isolated contributes to {}, short vertices contributes to {}".format(
    flag_explain_result[0],flag_explain_result[1],flag_explain_result[2])) # 826039, 199362, 3224
print("{} can not be explained by the three".format(sum(np.sum(unique_flag_array,axis=1)==0)))# 204994

# check how many miss belongs to basic missing set
ref_cause = list(set(problem_id_set).intersection(set(unique_gene_name_list)))
mut_cause = list(set(problem_id_set).difference(set(unique_gene_name_list)))
print("{} caused by ref and {} caused by mut".format(len(ref_cause),len(mut_cause)))

#
explainable_num = np.sum([np.sum(meta_dict[item]) > 0 for item in mut_cause])
print("{}/{} can be further explained".format(explainable_num,len(mut_cause)))
