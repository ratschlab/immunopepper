import numpy.random as npr

nucleotides_dict = {"C": '0', "G": '1', "A": '2', "T": '3'}
mut_options = {"0|0": '0', "0|1": '1', "1|0": '2', "1|1": '3'}

file = '../data_simulated/simulated_Ipp/variants.vcf'
new_file = '../data_simulated/data/variants.vcf'
#Read the vcf file

with open(file, 'r') as f:
    with open(new_file, 'w') as f_new:
        for line in f:
             if line.startswith('##'):
                 f_new.write(line)
                 continue
             if line.startswith('#'):
                fields = line.strip().split('\t')
                vcf_sample_set = fields[9:]
                #new_vcf = [name.split('/')[1].replace("_", "") for name in vcf_sample_set]
                new_vcf = vcf_sample_set
                fields[9:] = new_vcf
                line = '\t'.join(fields)
                f_new.write(line)
                f_new.write('\n')

                continue

             if line.startswith('1'):
                 line = line.strip().split('\t')
                 original = line[3]
                 index = nucleotides_dict.get(original)
                 random_idx = npr.randint(0, 4)
                 while random_idx == index:
                     random_idx = npr.randint(0, 4)

                # #Check the key to which the random_idx corresponds
                 for key, value in nucleotides_dict.items():
                     if value == str(random_idx):
                         line[4] = key

                 for i in range(9, len(line)):
                     random_index = npr.randint(0, 4)
                     for key, value in mut_options.items():
                            if value == str(random_index):
                                line[i] = key + ':0'

                 line_new = '\t'.join(line)

                # #Write the new line to the file overwriting the old line
                 prob = npr.random()
                 if prob < 0.005:
                    f_new.write(line_new)
                    f_new.write('\n')











