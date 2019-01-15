from immunopepper.constant import NOT_EXIST
import argparse
import sys
def parse_arguments(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", help="specify the directory that output data is in", required=False, default='')
    parser.add_argument("--mutation_mode", help="specify the mutation mdoe", required=False, default='ref')
    if len(argv) < 2:
        parser.print_help()
        sys.exit(1)
    pargs = parser.parse_args(argv)
    return pargs

def filter_kmer_dict_with_threshold(_dict,th=0):
    return {k:v for k,v in _dict.items() if v > th}

def build_kmer_dict(kmer_file):
    f = open(kmer_file,'r')
    kmer_dict = {}
    for i,line in enumerate(f):
        if i % 1000000 == 0:
            print("Processed {} lines".format(i))
        items = line.strip().split('\t')
        if len(items) == 2: # ignore abnormal case
            kmer = items[0]
            expr = float(items[1]) if NOT_EXIST != items[1] else NOT_EXIST
            if kmer not in kmer_dict:
                kmer_dict[kmer] = expr
            elif expr > kmer_dict[kmer]:
                kmer_dict[kmer] = expr
    return kmer_dict

def union_kmer_dict(_dict1,_dict2):
    union_key_set = set(_dict1).union(_dict2)
    new_dict = {}
    for key in union_key_set:
        if key in _dict1 and key in _dict2:
            expr = max(_dict1[key],_dict2[key])
        if key in _dict1 and key not in _dict2:
            expr = _dict1[key]
        if key not in _dict1 and key in _dict2:
            expr = _dict2[key]
        new_dict[key] = expr
    return new_dict

if __name__ == "__main__":
    arg = parse_arguments(sys.argv[1:])

    mutation_mode = arg.mutation_mode
    data_dir = arg.data_path

    back_kmer_file = data_dir+'{}_back_kmer.txt'.format(mutation_mode)
    junc_kmer_file = data_dir+'{}_junction_kmer.txt'.format(mutation_mode)
    concat_kmer_file = data_dir+'{}_concat_kmer.txt'.format(mutation_mode)

    # find kmer whose expression excess certain threshold
    back_kmer_dict = build_kmer_dict(back_kmer_file)
    junc_kmer_dict = build_kmer_dict(junc_kmer_file)
    concat_kmer_dict = build_kmer_dict(concat_kmer_file)

    full_kmer_dict = union_kmer_dict(junc_kmer_dict,concat_kmer_dict)
    filter_full_kmer = filter_kmer_dict_with_threshold(full_kmer_dict,0)
    neo_kmer_dict = set(filter_full_kmer).difference(back_kmer_dict)

    neo_kmer_file = data_dir+'{}_neo_kmer.txt'.format(mutation_mode)
    neo_kmer_file_fp = open(neo_kmer_file,'w')

    new_line_list = [neo_key + '\t'+str(full_kmer_dict[neo_key])+'\n' for neo_key in neo_kmer_dict]
    neo_kmer_file_fp.writelines(new_line_list)


