"""
Integrate multiple background kmer files and generate background kmer file
"""
def immunopepper_makebg(arg):
    kmer_file_list = arg.kmer_files_list
    output_file_path = arg.output_file_path
    verbose = arg.verbose
    uniq_kmer_set = set()
    for kmer_file in kmer_file_list:
        if verbose > 1:
            print(kmer_file)
        with open(kmer_file,'r') as f:
            kmer_list = [line.split('\t')[0] for line in f]
            uniq_kmer_set = uniq_kmer_set.union(kmer_list)
    uniq_kmer_list = sorted(uniq_kmer_set)
    with open(output_file_path,'w') as f_out:
        f_out.writelines("%s\n" % l for l in uniq_kmer_list)
    if verbose > 0:
        print("generate unique background kmer file in {}".format(output_file_path))


