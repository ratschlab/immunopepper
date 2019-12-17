from immunopepper.constant import NOT_EXIST
import argparse
import sys
import pickle
import os
def parse_arguments(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument("--back_kmer_dir", help="specify the directory that background kmer file is in", required=False, default='')
    parser.add_argument("--junc_kmer_dir", help="specify the directory that junction kmer file is in, also where the neokmer file will save in", required=False, default='')
    parser.add_argument("--concat_kmer_dir", help="specify the directory that concatenate kmer file is in", required=False, default='')
    parser.add_argument("--mutation_mode", help="specify the mutation mdoe", required=False, default='ref')
    if len(argv) < 2:
        parser.print_help()
        sys.exit(1)
    pargs = parser.parse_args(argv)
    return pargs

def filter_kmer_dict_with_threshold(_dict,th=0):
    return {k:v for k,v in list(_dict.items()) if v > th}

def build_kmer_dict(kmer_file):
    with open(kmer_file,'r') as f:
        kmer_dict = {}
        for i,line in enumerate(f):
            if i % 1000000 == 0:
                print("Processed {} lines".format(i))
            items = line.strip().split('\t')
            if len(items) == 3: # ignore abnormal case
                kmer = items[0]
                location = items[1]
                expr = float(items[2]) if NOT_EXIST != items[2] else NOT_EXIST
                if kmer not in kmer_dict:
                    kmer_dict[kmer] = [(location,expr)]
                else:
                    kmer_dict[kmer].append((location,expr))
    return kmer_dict

def union_kmer_dict(_dict1,_dict2):
    union_key_set = set(_dict1).union(_dict2)
    new_dict = {}
    for key in union_key_set:
        if key in _dict1 and key in _dict2:
            detail_list = _dict1[key]+_dict2[key]
        if key in _dict1 and key not in _dict2:
            detail_list = _dict1[key]
        if key not in _dict1 and key in _dict2:
            detail_list = _dict2[key]
        new_dict[key] = detail_list
    return new_dict

if __name__ == "__main__":
    arg = parse_arguments(sys.argv[1:])

    mutation_mode = arg.mutation_mode
    back_kmer_file = os.path.join(arg.back_kmer_dir,'{}_back_kmer.txt'.format(mutation_mode))
    junc_kmer_file = os.path.join(arg.junc_kmer_dir,'{}_junction_kmer.txt'.format(mutation_mode))
    concat_kmer_file = os.path.join(arg.concat_kmer_dir,'{}_concat_kmer.txt'.format(mutation_mode))

    # find kmer whose expression excess certain threshold
    back_kmer_dict = build_kmer_dict(back_kmer_file)
    junc_kmer_dict = build_kmer_dict(junc_kmer_file)
    concat_kmer_dict = build_kmer_dict(concat_kmer_file)

    full_kmer_dict = union_kmer_dict(junc_kmer_dict,concat_kmer_dict)
    filter_full_kmer = filter_kmer_dict_with_threshold(full_kmer_dict,0)
    filter_back_kmer = filter_kmer_dict_with_threshold(back_kmer_dict,0)
    neo_kmer_dict = set(filter_full_kmer).difference(filter_back_kmer)

    neo_kmer_file = os.path.join(arg.junc_kmer_dir,'{}_neo_kmer.txt'. format(mutation_mode))
    neo_kmer_file_fp = open(neo_kmer_file,'w')

    for neo_key in neo_kmer_dict:
        detail_list = full_kmer_dict[neo_key]
        detail_location = ';'.join([item[0] for item in detail_list])
        detail_expr = ';'.join([str(item[1]) for item in detail_list])
        new_line = '\t'.join([neo_key,detail_expr,detail_location])+'\n'
        neo_kmer_file_fp.write(new_line)


