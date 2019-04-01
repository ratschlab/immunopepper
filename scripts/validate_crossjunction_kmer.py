# validation work in March 9th
# The general purpose of this script is to generate the result
# published in cancer cell paper with ImmunoPepper
# Link: (https://www.sciencedirect.com/science/article/pii/S1535610818303064)

# validate 2019/03/07
import gzip
import pickle
import argparse
import sys
import os

import numpy as np
import pandas as pd
from immunopepper.constant import NOT_EXIST
from functools import reduce


def find_diff_coord(imm_gene_coord_dict,cancercell_gene_coord_dict):
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
        if gene_name not in cancercell_gene_coord_dict:
            cancercell_coord_set = {}
        else:
            cancercell_coord_set = cancercell_gene_coord_dict[gene_name]
        for coord,vertex_id in imm_coord_set:
            similar_result,found = fuzzy_comp_str2list(coord,cancercell_coord_set)
            if not found:
                gene_name_coord_str = gene_name+'_'+'_'.join([str(icoord) for icoord in coord])+'_'+vertex_id
                new_or_append_value_to_dict_key(unique_imm_coord_dict,gene_name,gene_name_coord_str)
    return unique_imm_coord_dict


def get_cancercell_ref_junction_dict(cancercell_junction_file):
    gt_lines = open(cancercell_junction_file,'r').readlines()
    ref_junction_dict = {} # (gene_name,coord_str_tuple) |-> peptide
    cancercell_gene_coord_dict = {} # gene_name |-> coord_str_tuple
    i = 0
    while i < len(gt_lines)-1:
        headline = gt_lines[i].strip()
        items = headline.split('_')
        gene_name = items[2]
        coord_str_tuple = tuple([int(coord)-1 if a%2 ==0 else int(coord) for a, coord in enumerate(items[-4:])]) # remove coord correction
        i += 1
        peptide = gt_lines[i].strip()
        ref_junction_dict[(gene_name,coord_str_tuple)] = peptide
        new_or_append_value_to_dict_key(cancercell_gene_coord_dict,gene_name,coord_str_tuple)
        i += 1
    return ref_junction_dict,cancercell_gene_coord_dict


def get_cancercell_kmer_dict(cancercell_result_file):
    def get_kmer_list(_str, k):
        if len(_str) < k:
            return [_str]
        else:
            kmer_list = [_str[i:i + k] for i in range(0, max(0, len(_str) - k + 1))]
            return kmer_list

    cancercell_pep_lines = open(cancercell_result_file, 'r').readlines()
    i = 0
    cancercell_kmer_dict = {}
    print("Parsing cancercell kmer result file. Need around 1 minute.")
    while i < len(cancercell_pep_lines):
        line = cancercell_pep_lines[i]
        gene_name = line.strip().split('_')[2]
        i += 1
        pep = cancercell_pep_lines[i].strip()
        kmer_list = get_kmer_list(pep, 9)
        new_dict = {kmer: gene_name for kmer in kmer_list}
        cancercell_kmer_dict = merge_two_dicts(cancercell_kmer_dict, new_dict)
        i += 1
    return cancercell_kmer_dict

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
    meta_df = pd.read_csv(meta_file,sep='\t')
    meta_flag_dict_key_with_coord = {}  # (geneName_fourCoords_vertexId) |-> flag_tuple
    meta_flag_dict_key_without_coord = {}
    imm_gene_coord_dict = {}  # gene_name |->List[(coord_str_tuple,vertex_id)], since one gene_name can have multuple lines
    for i in range(len(meta_df)):
        gene_name = meta_df['gene_name'][i]
        stop_flag = int(meta_df['has_stop_codon'][i])
        isolated_flag = int(meta_df['is_isolated'][i])
        som_variant_comb_num = int(len(meta_df['variant_comb'][i].split(';')) > 1)
        exon_coord = meta_df['exons_coor'][i]
        coord_str_tuple = tuple([int(coord) if coord != NOT_EXIST else NOT_EXIST for coord in exon_coord.split(';')])
        vertex_id = meta_df['vertex_idx'][i]
        key_with_coord = gene_name + '_' + '_'.join(exon_coord.split(';')) + '_' + vertex_id
        key_without_coord = gene_name+'_'+vertex_id.split(',')[0]+'_'+vertex_id.split(',')[1]
        start_v1 = coord_str_tuple[0]
        stop_v1 = coord_str_tuple[1]
        is_less_than3_flag = int((int(stop_v1) - int(start_v1)) < 3)
        flag_tuple = (stop_flag, isolated_flag, is_less_than3_flag,som_variant_comb_num)
        meta_flag_dict_key_with_coord = new_or_append_value_to_dict_key(meta_flag_dict_key_with_coord,key_with_coord,flag_tuple)
        meta_flag_dict_key_without_coord = new_or_append_value_to_dict_key(meta_flag_dict_key_without_coord,key_without_coord,flag_tuple)
        imm_gene_coord_dict = new_or_append_value_to_dict_key(imm_gene_coord_dict,gene_name,(coord_str_tuple, vertex_id))
    return meta_flag_dict_key_with_coord, meta_flag_dict_key_without_coord, imm_gene_coord_dict

def create_reference_cause_data():
    """
    Generate reference cause dictionary.
    It takes some time to generate so we will store those
    require dictionaries as pickle for subsequent use.

    Implement a comprehensive comparison between ImmunoPepper's reference_junction
    and CancerCell's reference_junction. We want to see how many outputs only
    exist in Immunopepper (denote as additional) and how many outputs only exist
    in CancerCell (denote as missing).
    """

    # load and parse immunopepper result
    imm_ref_meta_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/TCGA-13-1489/ref_metadata.tsv.gz'
    ref_meta_flag_dict_with_coord, ref_meta_flag_dict_without_coord, imm_gene_coord_dict = get_immunopepper_meta_dict(
        imm_ref_meta_file)

    # load and parse cancercell result
    cancercell_junction_file ='/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/' \
                               'analysis_pancan/ccell_rerun_2018/output/peptides.clean/split/REFERENCE.cj.fa'
    ref_junction_dict, cancercell_gene_coord_dict = get_cancercell_ref_junction_dict(cancercell_junction_file)

    # get the additional kmer dict (only appear in immunopepper result)
    additional_imm_coord_dict = find_diff_coord(imm_gene_coord_dict, cancercell_gene_coord_dict)

    # remove isolated case and only keep coord_str_tuple for further use
    imm_gene_coord_dict_without_vid = {
    gene_name: [coord_vid[0] for coord_vid in filter((lambda x: x[0][2] != NOT_EXIST), coord_vid_tuple_list)] for
    gene_name, coord_vid_tuple_list in list(imm_gene_coord_dict.items())}

    # get the missing kmer dict (only appear in cancercell paper result)
    missing_cancercell_gene_coord_dict = {
    gene_name: list(filter((lambda x: x in imm_gene_coord_dict_without_vid[gene_name]), coord_tuple_list)) for
    gene_name, coord_tuple_list in list(cancercell_gene_coord_dict.items())}

    # we know some of the missing kmers are due to extrapolation
    # which means there are 3 identical coords and only one coord is different.
    # We mark as True for those cases that can be explained by this.
    unique_cancercell_gene_explain_dict = {
    gene_name: [3 in np.sum((np.array(coord_tuple) - np.array(imm_gene_coord_dict_without_vid[gene_name])) == 0, axis=1)
                for coord_tuple in coord_tuple_list] for gene_name, coord_tuple_list in
    list(missing_cancercell_gene_coord_dict.items())}
    missing_explain_num = np.sum([sum(item) for item in list(unique_cancercell_gene_explain_dict.values())])  # 65360/65360

    # do a simple analysis of the result
    additional_junc_pair_num = np.sum([len(item) for item in list(additional_imm_coord_dict.values())])  # 1103270
    missing_junc_pair_num = np.sum([len(item) for item in list(missing_cancercell_gene_coord_dict.values())])  # 65360
    total_junc_pair_num = np.sum([len(item) for item in list(imm_gene_coord_dict.values())])  # 1622126
    total_cancercell_junc_pair_num = np.sum([len(item) for item in list(cancercell_gene_coord_dict.values())])  # 525539

    print("{} additional junc pair and {} miss junc pair in {} immuno total junc pairs"
          " and {} cancercell junc pairs. {} miss kmer are caused by extrapolation".format(additional_junc_pair_num,
                                                                                    missing_junc_pair_num,
                                                                                    total_junc_pair_num,
                                                                                    total_cancercell_junc_pair_num,
                                                                                    missing_explain_num))
    # try to explain the additional output
    additional_flag_list = []
    additional_gene_name_list = []
    for gene_name, additional_junc_coord_list in list(additional_imm_coord_dict.items()):
        additional_flag_list.extend(
            [reduce((lambda x, y: np.logical_or(x, y)), ref_meta_flag_dict_with_coord[additional_junc_coord]) for additional_junc_coord in
             additional_junc_coord_list])
        for additional_junc_coord in additional_junc_coord_list:
            items = additional_junc_coord.split('_')
            vertex_id = items[-1].split(',')
            gene_name = items[0]
            gene_name_str = '_'.join((gene_name, vertex_id[0], vertex_id[1]))
            additional_gene_name_list.append(gene_name_str)

    flag_explain_result = np.sum(np.array(additional_flag_list), axis=0)
    print("stop codon constributes to {}, isolated contributes to {}, short vertices contributes to {}".format(
        flag_explain_result[0], flag_explain_result[1], flag_explain_result[2]))  # 826039, 199362, 3224
    print("{} can not be explained by the three".format(sum(np.sum(np.array(additional_flag_list), axis=1) == 0)))  # 204994

    # load and parse reference junction kmer for further use
    # It can be used to explain some missing kmers.
    imm_ref_junction_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/TCGA-13-1489/ref_junction_kmer.txt'
    ref_kmer_dict = {line.split('\t')[0]:line.split('\t')[1] for line in open(imm_ref_junction_file,'r') if line.strip().split('\t')[-1] == 'True'}
    ref_kmer_set = set(ref_kmer_dict.keys())

    with open('ref_cause_dict2.pkl', 'wb') as f:
        pickle.dump((ref_kmer_set, imm_gene_coord_dict, cancercell_gene_coord_dict, additional_imm_coord_dict,
                      additional_gene_name_list, missing_cancercell_gene_coord_dict), f)


def parse_arguments(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", help="the sample names(can be string or ), can specify more than one sample", required=False, default='TCGA-13-1489')
    parser.add_argument("--mutation_mode", help="specify the mutation mdoe", required=False, default='germline')
    if len(argv) < 2:
        parser.print_help()
        sys.exit(1)

    pargs = parser.parse_args(argv)
    return pargs


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
cancercell_kmer_file = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/analysis_pancan/' \
                   'ccell_rerun_2018/output/peptides.clean/split/cj_kmers/{}.cj.{}.cj_kmers_9.fa'.format(sample_name,mutation_mode)
mutation_meta_gz_file = '/cluster/work/grlab/projects/TCGA/immunopepper_rerun/{}/{}_metadata.tsv.gz'.format(sample_name,mutation_mode)

if mutation_mode == 'somatic_and_germline':
    cancercell_kmer_file = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/analysis_pancan/' \
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
if not os.path.exists('ref_cause_dict.pkl'):
    print("Reference cause dictionary does not exist. Begin generating and might take 20 minutes...\n")
    create_reference_cause_data()
f = open('ref_cause_dict.pkl','rb')
ref_kmer_set,imm_gene_coord_dict,cancercell_gene_coord_dict,additional_imm_coord_dict,additional_gene_name_list,missing_cancercell_gene_coord_dict = pickle.load(f)

# get immunopepper junction metadata dict
_,mut_meta_flag_dict_without_coord,_ = get_immunopepper_meta_dict(mutation_meta_gz_file)
mut_meta_flag_dict_without_coord = merge_two_dicts(mut_meta_flag_dict_without_coord,aux_mut_meta_flag_dict_without_coord)

# get cancercell kmer dict
# kmer: -> gene_name
cancercell_kmer_dict = get_cancercell_kmer_dict(cancercell_kmer_file)
cancercell_kmer_set = set(cancercell_kmer_dict.keys())

# get immunopepper kmer dict, only focus on the cross junction part
# kmer: -> gene_name
immunopepper_dict = {line.split('\t')[0]:line.split('\t')[1] for line in open(immunopepper_file,'r') if line.strip().split('\t')[-1] == 'True'}
immunopepper_dict = merge_two_dicts(immunopepper_dict,aux_germ_immunopepper_dict)

########
# Part 1: get kmer dict
########

# get all kmer returned by immunopepper
immunopepper_kmer = set(immunopepper_dict.keys())
common_kmer_set = cancercell_kmer_set.intersection(immunopepper_kmer)
additional_kmer_list = list(immunopepper_kmer.difference(cancercell_kmer_set))
missing_kmer_list = list(cancercell_kmer_set.difference(immunopepper_kmer))
num_common_kmer = len(common_kmer_set)
num_additional_kmer = len(additional_kmer_list)
num_missing_kmer = len(missing_kmer_list)
s_summary = ">>>>>>>>Validation Start\n\nComparison overview. (additional kmer means imm subtracts mat, miss kmer means mat subtracts imm).\n" \
            ">> {} common kmers, {} additional kmers and {} miss kmers".format(num_common_kmer,num_additional_kmer,num_missing_kmer)
print(s_summary)

########
# Part 2: Explain missing kmer
########

# case 1: Mutations occur at the extrapolation part
#   so the whole peptide is kept alghough they are the same with ref peptide near the cross junction place
#   these missing kmers are not important because they will be (hopefully) covered by the ref_kmer_set
# case 2: Extrapolation put forward the promotor and thus might generate more kmers. These missing kmer
#   actually come from the missing_cancercell_gene_coord_dict.

# consider case 1
unexplained_missing_kmer_list = list(set(missing_kmer_list).difference(ref_kmer_set))
# consider case 3
num_can_not_explained_missing_kmer = np.sum([cancercell_kmer_dict[missing_kmer] not in missing_cancercell_gene_coord_dict for missing_kmer in unexplained_missing_kmer_list])

s_explain_missing = "Explain the {} miss kmers, ideally they can all be found in reference kmer list.\n" \
                    ">> {} kmers can not be found in reference kmer".format(num_missing_kmer,num_can_not_explained_missing_kmer)
print(s_explain_missing)


########
# Part 2: Explain additional kmer
########

# get the corresponding gene_name str for those additional kmer
# we expect most parts of it should also appear in the additional junction kmer dict
# that is additional_gene_name_list.
additional_kmer_id_list = list(set([immunopepper_dict[_kmer] for _kmer in additional_kmer_list]))

########
# Part 2.1: Consider the reference cause first
#   we compare cancercell's reference junction output with immunopepper's,
#   and see how many differences we have.
#   Most of the difference in kmer comparison can attribute to difference in junction
########

ref_cause = list(set(additional_kmer_id_list).intersection(set(additional_gene_name_list)))
mut_cause = list(set(additional_kmer_id_list).difference(set(additional_gene_name_list)))
s_explain_addition = "Explain the {} additional kmers. Remember one junction output can generate" \
                     " 9 kmers. We can see how many difference " \
                     "attribute to reference kmer and how many to mutatation-specific kmer.\n" \
                     ">> {} caused by ref and {} caused by mut in the total {} " \
                     "additional junction output".format(num_additional_kmer,len(ref_cause),len(mut_cause),len(additional_kmer_id_list))
print(s_explain_addition)

########
# Part 2.2: Consider the mutation cause
########
# case 1: junction pair appears in both mat and imm. However, germline/somatic mutation introduces a stop codon and
#   mat will not output the peptide but imm will do.
# case 2: Multiple somatic may introduce changes over the cross_junction region and thus create new kmers
# case 3: complicated case. ENSG00000034152.14_34_35
#   In this gene, we have the path 24-26-30-34-35/38
#   germline mutation introduces one stop codon in (24,26), so mat missed all the following pairs
#   (26,30),(30,34),(34,35),(34,38).
#   but we do not see (26,30) or (30,34) in immunopeppper_dict (thus not see in mut_case) why?
#   because (26,31) junction pair has exactly the same kmer output as (26,30)
#   since we will cover the original mapping of (kmer->(26,30)) and change it to be (kmer->(26,31))
#   So ony (26,31) are left in additional_kmer_id_list.
#   since (26,31) has stop codon so they can be explained by ref_cause and not seen in mut_case.
#   Therefore, we finally see (24,26), (34,35), (34,38) in mut_case
#   (24,26) has stop codon and can be further explained but (34,35) and (34,38) is the difference
#   caused by (24,26) and have no flag so they can not be explicitly explained.
#   But we can see, this case is very rare for single sample+mutation type.

explainable_num = np.sum([np.sum(mut_meta_flag_dict_without_coord[item]) > 0 for item in mut_cause])
s_explain_mut_cause = "Explain the mut causes:\n" \
                      ">> {}/{} can be further explained".format(explainable_num,len(mut_cause))
print(s_explain_mut_cause)

########
# Part 3: Summary
########
# let's see how these different reasons contribute to the additional kmer
additional_flag_tuple = np.array([reduce((lambda x,y: np.logical_or(x,y)),mut_meta_flag_dict_without_coord[item]) for item in additional_kmer_id_list])
flag_explain = np.sum(additional_flag_tuple,axis=0)

# The effect of extrapolation is not so easy to get directly
# Here we just do a subtraction.
num_extrapolation_explain = np.sum(np.sum(additional_flag_tuple,axis=1) == 0)
num_can_not_explained_additional_kmer = len(mut_cause)-explainable_num

s_final_conclusion = "\n>>>>>>>>Valiadtion Summary for {} {}\nThere are {} common kmers, {} missing kmers only appear " \
                     "in Matthias's result" \
                     ", {} additional kmers only appear in Immunopepper's result. {} missing kmers and {} additional kmers can not be easily " \
                     "explained.\n For the {} concerned additional junction, {} are caused by stop codon, " \
                     "{} by isolated, {} by short vertices, {} by somatic variant combination, {} by extrapolation.".format(
                    sample_name,mutation_mode,num_common_kmer,num_missing_kmer,num_additional_kmer,
                                        num_can_not_explained_missing_kmer,num_can_not_explained_additional_kmer,
                    len(additional_kmer_id_list), flag_explain[0], flag_explain[1], flag_explain[2],flag_explain[3],num_extrapolation_explain)
print(s_final_conclusion)


# The following code are kept for further convenience
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
