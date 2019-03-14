# validation work in March 9th
# validate 03/07/2019
import numpy as np
import gzip
import cPickle
import argparse
import sys

from immunopepper.constant import NOT_EXIST

def get_immunopepper_meta_dict(meta_file):
    """
    Why we have with_coord and without_coord?
    It's due to the current mismatch. Immunopepper peptide.fa only has gene_name+vertex_id
    while mat's peptide.fa has gene_name+vertex.
    To try to explain the difference of reference peptide output,
    we need build with_coord_dict.
    However, when we trace back the flag tuple of those additional kmer to explain them,
    we need without_coord dict.
    In summary, to generate reference comparison result, we need with_coord dict;
    for validating other samples and mutation types, we only need without_coord_dict
    """
    f_meta = gzip.open(meta_file, 'r')
    f_meta.readline()
    meta_flag_dict_key_with_coord = {}  # (geneName_fourCoords_vertexId) |-> flag_tuple
    meta_flag_dict_key_without_coord = {}
    imm_gene_coord_dict = {}  # gene_name |->List[(coord_str_tuple,vertex_id)], since one gene_name can have multuple lines
    for line in f_meta:
        items = line.strip().split('\t')
        gene_name = items[2]
        stop_flag = int(items[9])
        isolated_flag = int(items[11])
        som_variant_comb_num = int(len(items[12].split(';')) > 1)
        coord_str_tuple = tuple([int(coord) if coord != NOT_EXIST else NOT_EXIST for coord in items[14].split(';')])
        vertex_id = items[15]
        key_with_coord = gene_name + '_' + '_'.join(items[14].split(';')) + '_' + vertex_id
        key_without_coord = gene_name+'_'+vertex_id.split(',')[0]+'_'+vertex_id.split(',')[1]
        start_v1 = coord_str_tuple[0]
        stop_v1 = coord_str_tuple[1]
        is_less_than3_flag = int((int(stop_v1) - int(start_v1)) < 3)
        flag_tuple = (stop_flag, isolated_flag, is_less_than3_flag,som_variant_comb_num)
        meta_flag_dict_key_with_coord = new_or_append_value_to_dict_key(meta_flag_dict_key_with_coord,key_with_coord,flag_tuple)
        meta_flag_dict_key_without_coord = new_or_append_value_to_dict_key(meta_flag_dict_key_without_coord,key_without_coord,flag_tuple)
        imm_gene_coord_dict = new_or_append_value_to_dict_key(imm_gene_coord_dict,gene_name,(coord_str_tuple, vertex_id))
    return meta_flag_dict_key_with_coord, meta_flag_dict_key_without_coord, imm_gene_coord_dict


def parse_arguments(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", help="the sample names(can be string or ), can specify more than one sample", required=False, default='TCGA-13-1489')
    parser.add_argument("--mutation_mode", help="specify the mutation mdoe", required=False, default='germline')
    if len(argv) < 2:
        parser.print_help()
        sys.exit(1)

    pargs = parser.parse_args(argv)
    return pargs


def get_kmer_list(_str,k):
    if len(_str) < k:
        return [_str]
    else:
        kmer_list = [_str[i:i+k] for i in range(0,max(0,len(_str)-k+1))]
        return kmer_list


def new_or_append_value_to_dict_key(_dict,key,value):
    if key in _dict:
        _dict[key].append(value)
    else:
        _dict[key] = [value]
    return _dict


def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

arg = parse_arguments(sys.argv[1:])
sample_name = arg.samples
mutation_mode = arg.mutation_mode
immunopepper_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/{}/{}_junction_kmer.txt'.format(sample_name,mutation_mode)
nora_result_file = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/analysis_pancan/' \
                   'ccell_rerun_2018/output/peptides.clean/split/cj_kmers/{}.cj.{}.cj_kmers_9.fa'.format(sample_name,mutation_mode)
mutation_meta_gz_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/{}/{}_metadata.tsv.gz'.format(sample_name,mutation_mode)

if mutation_mode == 'somatic_and_germline':
    nora_result_file = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/analysis_pancan/' \
                       'ccell_rerun_2018/output/peptides.clean/split/cj_kmers/{}.cj.{}.cj_kmers_9.fa'.format(
        sample_name, 'germline_somatic')
    aux_germ_immunopepper_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/{}/{}_junction_kmer.txt'.format(
        sample_name, 'germline')
    aux_germ_immunopepper_dict = {line.split('\t')[0]:line.split('\t')[1] for line in open(aux_germ_immunopepper_file,'r') if line.strip().split('\t')[-1] == 'True'}
    aux_germ_meta_gz_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/{}/{}_metadata.tsv.gz'.format(sample_name,'germline')
    _, aux_mut_meta_flag_dict_without_coord, _ = get_immunopepper_meta_dict(aux_germ_meta_gz_file)

else:
    # In Immunopepper's implementation, in the somatic_and_germline mode (in germline_somatic Immunopepper),
    # the reference peptide is the germline-applied while in Matthias', the reference peptide is just the reference.
    # Therefore, quite a lot of kmers are already included in 'germline_junction_kmer' and not included in
    # 'somatic_germline_junction_kmer'. In Matthias' output, those germline kmers are also included in 'germline_somatic_cj_kmers'
    # To compare the two, we need the auxillary dict generated from germline mode.
    aux_germ_immunopepper_dict = {}
    aux_mut_meta_flag_dict_without_coord = {}

########
# Part 0:  load some auxillary function and data
########
'''
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
ENSG00000076344.11_12_11
# load sequence file
import Bio.SeqIO as BioIO
ref_path = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/sequence/genome.fa'
seq_dict = {}
interesting_chr = map(str, range(1, 23)) + ["X", "Y", "MT"]
print('Parsing genome sequence ...')
for record in BioIO.parse(ref_path, "fasta"):
    if record.id in interesting_chr:
        seq_dict[record.id] = str(record.seq).strip()
'''

f = open('ref_cause_dict.pkl','rb')
ref_kmer_set,imm_gene_coord_dict,mat_gene_coord_dict,unique_imm_coord_dict,unique_gene_name_list,unique_mat_gene_coord_dict = cPickle.load(f)

_,mut_meta_flag_dict_without_coord,_ = get_immunopepper_meta_dict(mutation_meta_gz_file)
mut_meta_flag_dict_without_coord = merge_two_dicts(mut_meta_flag_dict_without_coord,aux_mut_meta_flag_dict_without_coord)

# parse nora's result from 17mer to 9mer
nora_pep_lines = open(nora_result_file,'r').readlines()
i = 0
pep_dict = {}
nora_kmer_list = []
nora_kmer_dict = {}
while i < len(nora_pep_lines):
    if i % 10000 == 0:
        print(i)
    line = nora_pep_lines[i]
    gene_name = line.strip().split('_')[2]
    i += 1
    pep = nora_pep_lines[i].strip()
    kmer_list = get_kmer_list(pep,9)
    new_dict = {kmer:gene_name for kmer in kmer_list}
    nora_kmer_dict = merge_two_dicts(nora_kmer_dict,new_dict)
    nora_kmer_list.extend(kmer_list)
    i += 1
nora_kmer_set = set(nora_kmer_list)

##########
# Part 1: compare kmer output
##########
# construct immuno_pepper dict, only focus on the cross junction part
# kmer: -> immuno_pepper name
immunopepper_dict = {line.split('\t')[0]:line.split('\t')[1] for line in open(immunopepper_file,'r') if line.strip().split('\t')[-1] == 'True'}
immunopepper_dict = merge_two_dicts(immunopepper_dict,aux_germ_immunopepper_dict)

# get all kmer returned by immunopepper
immunopepper_kmer = set(immunopepper_dict.keys())
common_kmer_set = nora_kmer_set.intersection(immunopepper_kmer)
additional_kmer_list = list(immunopepper_kmer.difference(nora_kmer_set))
miss_kmer_list = list(nora_kmer_set.difference(immunopepper_kmer))
num_common_kmer = len(common_kmer_set)
num_additional_kmer = len(additional_kmer_list)
num_miss_kmer = len(miss_kmer_list)
s_summary = ">>>>>>>>Validation Start\n\nComparison overview. (additional kmer means imm subtracts mat, miss kmer means mat subtracts imm).\n" \
            ">> {} common kmers, {} additional kmers and {} miss kmers".format(num_common_kmer,num_additional_kmer,num_miss_kmer)
print(s_summary)

# check missing kmer
# case 1: Mutations occur at the extrapolation part
#   so the whole peptide is kept alghough they are the same with ref peptide near the cross junction place
#   these missing kmers are not important because they will be (hopefully) covered by the ref_kmer_set
# case 2: Extrapolation put forward the promotor and thus might generate more kmers. These missing kmer
#   actually come from the unique_mat_gene_coord_dict.
'''
immunopepper_ref_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/TCGA-13-1489/ref_junction_kmer.txt'
ref_kmer_dict = {line.split('\t')[0]:line.split('\t')[1] for line in open(immunopepper_ref_file,'r') if line.strip().split('\t')[-1] == 'True'}
ref_kmer_set = set(ref_kmer_dict.keys())
'''
unexplained_miss_kmer_list = list(set(miss_kmer_list).difference(ref_kmer_set))
num_can_not_explained_miss_kmer = np.sum([nora_kmer_dict[miss_kmer] not in unique_mat_gene_coord_dict for miss_kmer in unexplained_miss_kmer_list])

s_explain_missing = "Explain the {} miss kmers, ideally they can all be found in reference kmer list.\n" \
                    ">> {} kmers can not be found in reference kmer".format(num_miss_kmer,num_can_not_explained_miss_kmer)
print(s_explain_missing)

problem_id_list = list(set([immunopepper_dict[_kmer] for _kmer in additional_kmer_list]))

'''
# Some examples
from immunopepper.utils import translate_dna_to_peptide
from immunopepper.utils import complementary_seq
# example 1
# Extrapolation introduce stop codon
seq7 = seq_dict['7']
print(translate_dna_to_peptide(complementary_seq(seq7[2854510:2854891][::-1]))) # ('LSSLLDCRSADWLPP', True) Matthias implementation (start from extrapolation)
print(translate_dna_to_peptide(complementary_seq(seq7[2854510:2854564][::-1]))) # ('MKRRMFPRPCLARMPGSR', False) immunopepper implementation (start from cds)

# example 2 (one of the special 58)
# Mutation introduce stop codon
seq = seq_dict['8']
print(translate_dna_to_peptide(complementary_seq(seq[101724589:101724685][::-1]+seq[101721686:101721947][::-1]))) # ('LSSLLDCRSADWLPP', True) Matthias implementation (start from extrapolation)
print(translate_dna_to_peptide(complementary_seq(seq[101724589:101724685][::-1]))) # ('LSSLLDCRSADWLPP', True) Matthias implementation (start from extrapolation)
'''

#########
# Part 2: calculate the reference cause
# we compare matthias's reference junction output with immunopepper's,
# and see how many differences we have.
# Most of the difference in kmer comparison can attribute to difference in junction
#########
# build the missing table via reference comparing
'''
f_gt = open('/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/analysis_pancan/ccell_rerun_2018/output/peptides.clean/split/REFERENCE.cj.fa','r')
def get_matthias_ref_junction_dict(f_gt):
    gt_lines = f_gt.readlines()
    ref_junction_dict = {} # (gene_name,coord_str_tuple) |-> peptide
    mat_gene_coord_dict = {} # gene_name |-> coord_str_tuple
    i = 0
    while i < len(gt_lines)-1:
        headline = gt_lines[i].strip()
        items = headline.split('_')
        gene_name = items[2]
        coord_str_tuple = tuple([int(coord)-1 if a%2 ==0 else int(coord) for a, coord in enumerate(items[-4:])]) # remove coord correction
        i += 1
        peptide = gt_lines[i].strip()
        ref_junction_dict[(gene_name,coord_str_tuple)] = peptide
        new_or_append_value_to_dict_key(mat_gene_coord_dict,gene_name,coord_str_tuple)
        i += 1
    return ref_junction_dict,mat_gene_coord_dict

ref_meta_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/TCGA-13-1489/ref_metadata.tsv.gz'

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
                new_or_append_value_to_dict_key(unique_imm_coord_dict,gene_name,gene_name_coord_str)
    return unique_imm_coord_dict

ref_meta_flag_dict_with_coord,ref_meta_flag_dict_without_coord, imm_gene_coord_dict = get_immunopepper_meta_dict(ref_meta_file)
ref_junction_dict, mat_gene_coord_dict = get_matthias_ref_junction_dict(f_gt)
unique_imm_coord_dict = find_diff_coord(imm_gene_coord_dict,mat_gene_coord_dict)

# explore the missing junction
# ideally all of them are caused by extrapolation
imm_gene_coord_dict_without_vid = {gene_name:[coord_vid[0] for coord_vid in filter((lambda x: x[0][2] != NOT_EXIST),coord_vid_tuple_list)]for gene_name, coord_vid_tuple_list in imm_gene_coord_dict.items()}
unique_mat_gene_coord_dict = {gene_name:filter((lambda x: x in imm_gene_coord_dict_without_vid[gene_name]),coord_tuple_list) for gene_name,coord_tuple_list in mat_gene_coord_dict.items()}
unique_mat_gene_explain_dict = {gene_name:[3 in np.sum((np.array(coord_tuple)-np.array(imm_gene_coord_dict_without_vid[gene_name])) == 0,axis=1) for coord_tuple in coord_tuple_list] for gene_name,coord_tuple_list in unique_mat_gene_coord_dict.items()}
miss_explain_num = np.sum([sum(item) for item in unique_mat_gene_explain_dict.values()]) # 65360/65360

# do a simple analysis of the result
addition_junc_pair_num = np.sum([len(item) for item in unique_imm_coord_dict.values()]) # 1103270
miss_junc_pair_num = np.sum([len(item) for item in unique_mat_gene_coord_dict.values()]) # 65360
total_junc_pair_num = np.sum([len(item) for item in imm_gene_coord_dict.values()]) # 1622126
total_mat_junc_pair_num = np.sum([len(item) for item in mat_gene_coord_dict.values()]) # 525539

print("{} additional junc pair and {} miss junc pair in {} immuno total junc pairs"
      " and {} mat junc pairs. {} miss kmer are caused by extrapolation".format(addition_junc_pair_num,miss_junc_pair_num,total_junc_pair_num,total_mat_junc_pair_num,miss_explain_num))

# see detail information
unique_flag_list = []
#unique_gene_name_list = []
for gene_name, miss_junc_list in unique_imm_coord_dict.items():
    unique_flag_list.extend([reduce((lambda x,y: np.logical_or(x,y)),ref_meta_flag_dict_with_coord[miss_junc]) for miss_junc in miss_junc_list])
    for miss_junc in miss_junc_list:
        items = miss_junc.split('_')
        vertex_id = items[-1].split(',')
        gene_name = items[0]
        gene_name_str = '_'.join((gene_name,vertex_id[0],vertex_id[1]))
        unique_gene_name_list.append(gene_name_str)

with open('ref_cause_dict2.pkl','wb') as f:
    cPickle.dump((ref_kmer_set,imm_gene_coord_dict,mat_gene_coord_dict,unique_imm_coord_dict,unique_gene_name_list,unique_mat_gene_coord_dict),f)
unique_flag_array = np.array(unique_flag_list)
flag_explain_result = np.sum(unique_flag_array, axis=0)
print("stop codon constributes to {}, isolated contributes to {}, short vertices contributes to {}".format(
    flag_explain_result[0],flag_explain_result[1],flag_explain_result[2])) # 826039, 199362, 3224
print("{} can not be explained by the three".format(sum(np.sum(unique_flag_array,axis=1)==0)))# 204994
'''


# check how many miss belongs to basic missing set
ref_cause = list(set(problem_id_list).intersection(set(unique_gene_name_list)))
mut_cause = list(set(problem_id_list).difference(set(unique_gene_name_list)))
s_explain_addition = "Explain the {} additional kmers. Remember one junction output can generate" \
                     " 9 kmers. We can see how many difference " \
                     "attribute to reference kmer and how many to mutatation-specific kmer.\n" \
                     ">> {} caused by ref and {} caused by mut in the total {} " \
                     "additional junction output".format(num_additional_kmer,len(ref_cause),len(mut_cause),len(problem_id_list))
print(s_explain_addition)

#
# case 1: junction pair appears in both mat and imm. However, germline/somatic mutation introduces a stop codon and
#   mat will not output the peptide but imm will do.
# case 2: Multiple somatic may introduce changes over the cross_junction region and thus create new kmers
# case 3: complated case. ENSG00000034152.14_34_35
#   In this gene, we have the path 24-26-30-34-35/38
#   germline mutation introduces one stop codon in (24,26), so mat missed all the following pairs
#   (26,30),(30,34),(34,35),(34,38).
#   but we do not see (26,30) or (30,34) in immunopeppper_dict (thus not see in mut_case) why?
#   because (26,31) junction pair has exactly the same kmer output as (26,30)
#   since we will cover the original mapping of (kmer->(26,30)) and change it to be (kmer->(26,31))
#   So ony (26,31) are left in problem_id_list.
#   since (26,31) has stop codon so they can be explained by ref_cause and not seen in mut_case.
#   Therefore, we finally see (24,26), (34,35), (34,38) in mut_case
#   (24,26) has stop codon and can be further explained but (34,35) and (34,38) is the difference
#   caused by (24,26) and have no flag so they can not be explicitly explained.
#   But we can see, this case is very rare for single sample+mutation type.

explainable_num = np.sum([np.sum(mut_meta_flag_dict_without_coord[item]) > 0 for item in mut_cause])
addition_flag_tuple = np.array([reduce((lambda x,y: np.logical_or(x,y)),mut_meta_flag_dict_without_coord[item]) for item in problem_id_list])
flag_explain = np.sum(addition_flag_tuple,axis=0)
num_extrapolation_explain = np.sum(np.sum(addition_flag_tuple,axis=1) == 0)
s_explain_mut_cause = "Explain the mut causes:\n" \
                      ">> {}/{} can be further explained".format(explainable_num,len(mut_cause))
print(s_explain_mut_cause)
num_can_not_explained_additional_kmer = len(mut_cause)-explainable_num

s_final_conclusion = "\n\n>>>>>>>>>>>>>>>>>>Valiadtion Summary for {} {}\nThere are {} common kmers, {} missing kmers only appear " \
                     "in Matthias's result" \
                     ", {} additional kmers only appear in Immunopepper's result. {} missing kmers and {} can not be easily " \
                     "explained.\n For the {} concerned additional junction, {} are caused by stop codon, " \
                     "{} by isolated, {} by short vertices, {} by somatic variant combination, {} by extrapolation.".format(
                    sample_name,mutation_mode,num_common_kmer,num_miss_kmer,num_additional_kmer,
                                        num_can_not_explained_miss_kmer,num_can_not_explained_additional_kmer,
                    len(problem_id_list), flag_explain[0], flag_explain[1], flag_explain[2],flag_explain[3],num_extrapolation_explain)
print(s_final_conclusion)
