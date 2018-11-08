import os
import gzip
pos_result_path = '/Users/jiayu/PycharmProjects/CBM_RA/projects2018_immunopepper/tests/test1/test1neg'
neg_esult_path = '/Users/jiayu/PycharmProjects/CBM_RA/projects2018_immunopepper/tests/test1/test1pos'

def rename_as_gt_file(result_path):
    file_list = os.listdir(result_path)
    for file in file_list:
        item = file.split('.')
        suffix = item[-1]
        if suffix == 'gz':
            file_path = os.path.join(result_path, file)
            lines = gzip.open(file_path,'rb')
            new_file_name = item[0]+'_gt.tsv'
            new_file_path = os.path.join(result_path,new_file_name)
            with open(new_file_path, 'w') as ftsv:
                ftsv.writelines(lines.read())
        else:  # fa file
            name = item[0]
            name_part = name.split('_')
            if name_part[-1] != 'gt':
                name_part.append('gt')
                new_name = '_'.join(name_part)+'.'+suffix
                file_path = os.path.join(result_path,file)
                new_file_path = os.path.join(result_path,new_name)
                os.rename(file_path,new_file_path)
rename_as_gt_file(pos_result_path)
rename_as_gt_file(neg_esult_path)

# split the groundtruth to positive and negative part
# def change_neg_gene_id(neg_lines,mode='tsv'):
#     new_neg_lines = []
#     if mode == 'tsv':
#         for line in neg_lines:
#             items = line.split('\t')
#             id_items = items[0].split('.')
#             id_items[0] = '0'
#             new_id = '.'.join(id_items)
#             items[0] = new_id
#             new_neg_lines.append('\t'.join(items))
#     if mode == 'fa':
#         for i, line in enumerate(neg_lines):
#             if i%2 != 0:
#                 line_list = list(line)
#                 line_list[1] = '0'
#                 new_line = ''.join(line_list)
#                 new_neg_lines.append(new_line)
#             else:
#                 new_neg_lines.append(line)
#     return new_neg_lines
#
# for file in file_list:
#     item = file.split('.')
#     suffix = item[1]
#     name = item[0]
#     new_posname = name+'_pos'+'.'+suffix
#     new_negname = name+'_neg'+'.'+suffix
#     pos_path = os.path.join(result_path,new_posname)
#     neg_path = os.path.join(result_path,new_negname)
#     ori_path = os.path.join(result_path,file)
#     with open(pos_path,'w') as fpos, open(neg_path,'w') as fneg, open(ori_path,'r') as fori:
#         lines = fori.readlines()
#         if suffix == 'tsv':
#             assert ((len(lines) - 1) % 2 == 0)
#             mid = (len(lines) + 1) / 2
#             new_poslines = lines[0]+''.join(lines[1:mid])
#             neg_lines = change_neg_gene_id(lines[mid:],'tsv')
#             new_neglines = lines[0]+''.join(neg_lines)
#             fpos.writelines(new_poslines)
#             fneg.writelines(new_neglines)
#         else:
#             assert(len(lines)%2 == 0)
#             mid = len(lines)/2
#             new_poslines = lines[:mid]
#             new_neglines = change_neg_gene_id(lines[mid:], 'fa')
#             fpos.writelines(new_poslines)
#             fneg.writelines(new_neglines)
#     fpos.close()
#     fneg.close()




