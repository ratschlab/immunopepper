import argparse
import gzip
import re

def get_exon_dict(meta_file_path):
    if meta_file_path.endswith('gz'):
        f = gzip.open(meta_file_path, 'r')
    else:
        f = open(meta_file_path, 'r')
    exon_dict = {}
    for line in f:
        ### handle header
        if line.startswith('output_id'):
            header = line.strip().split('\t')
            eidx = header.index('exons_coor')
            continue
        items = line.strip().split('\t')
        exon_list = items[eidx].split(';')
        idx = items[0]
        gene_name = items[2]
        strand = items[4]
        if strand == '+':
            key = (gene_name, exon_list[1],exon_list[2])
            if key in exon_dict:
                exon_dict[key].append((idx,exon_list[0],exon_list[3]))
            else:
                exon_dict[key] = [(idx,exon_list[0],exon_list[3])]
        else:
            key = (gene_name, exon_list[0], exon_list[3])
            if key in exon_dict:
                exon_dict[key].append((idx, exon_list[2], exon_list[1]))
            else:
                exon_dict[key] = [(idx, exon_list[2], exon_list[1])]
    return exon_dict

def get_remove_id(exon_dict):
    remove_id_list = []
    for exon_pair_list in exon_dict.values():
        L = len(exon_pair_list)
        if L < 2:
            continue
        for i, exon_pair in enumerate(exon_pair_list):
            for j in range(L):
                i_pos1 = exon_pair[1]
                i_pos2 = exon_pair[2]
                j_pos1 = exon_pair_list[j][1]
                j_pos2 = exon_pair_list[j][2]
                if j!=i and i_pos1 >= j_pos1 and i_pos2 <= j_pos2:
                    if exon_pair[0] == '33803.7':
                        print exon_pair_list
                    remove_id_list.append(exon_pair[0])
                    break
    return remove_id_list


def get_filtered_output(meta_file, peptide_file, remove_id_list):
    if meta_file.endswith('gz'):
        fin1 = gzip.open(meta_file, 'r')
        fout1 = gzip.open(re.sub(r'.tsv.gz$', '', meta_file) + '.filt.tsv.gz', 'w')
    else:
        fin1 = open(meta_file, 'r')
        fout1 = open(re.sub(r'.tsv$', '', meta_file) + '.filt.tsv', 'w')
    fin2 = open(peptide_file)
    fout2 = open(re.sub(r'.fa$', '', peptide_file) + '.filt.fa', 'w')
    header_line = fin1.readline()
    fout1.write(header_line)
    lines1 = fin1.readlines()
    lines2 = fin2.readlines()
    for i,line in enumerate(lines1):
        items = line.strip().split('\t')
        if items[0] not in remove_id_list:
            print(items[0])
            fout1.write(line)
            fout2.write(lines2[2*i])
            fout2.write(lines2[2*i+1])
    fin1.close()
    fin2.close()
    fout1.close()
    fout2.close()

if "__name__" is not "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--meta_file", nargs=1, help="the generated meta tsv file name", required=True)
    parser.add_argument("--peptide_file", nargs=1, help="the corresponding peptide file", required=True)
    args = parser.parse_args()
    meta_file = args.meta_file[0]
    peptide_file = args.peptide_file[0]
    # meta_file = 't/TCGA-13-1489/ref_meta1.tsv.gz'
    # peptide_file = 't/TCGA-13-1489/ref_peptide2.fa'
    exon_dict = get_exon_dict(meta_file)
    remove_id_list = get_remove_id(exon_dict)
    get_filtered_output(meta_file, peptide_file, remove_id_list)




                    

