import numpy as np
import numpy.random as npr
import os
import gzip
path_junc = "../../../immunopepper_usecase/cohort_mutNone/ref_graph_kmer_JuncExpr"
path_normals = "../../../immunopepper_usecase/filter_case/ref_graph_kmer_NormalExpr"
path_unified_cancer = "../../../immunopepper_usecase/filter_case/ref_graph_kmer_CancerExpr"

#Fix npr seed
npr.seed(42)

#Create folders if they don't exist
if os.path.exists(path_normals) == False:
    os.mkdir(path_normals)
if os.path.exists(path_unified_cancer) == False:
    os.mkdir(path_unified_cancer)
first = 0
#We first drop some of the kmers and unify all the junction files into one file
# Read the Junction files in the path_junc folder
lines_cancer = []
lines_normal = []
for filename in os.listdir(path_junc):
    file_path = os.path.join(path_junc, filename)

    if os.path.isfile(file_path):
        with gzip.open(file_path, 'rt') as f:
            lines = f.readlines()
            kmer_info = lines[1:]
            if first == 0:
                lines_normal.append(lines[0])
                lines_cancer.append(lines[0])
                first = 1
            lines_cancer.append(lines[1:])

            for line in kmer_info:
                prob_drop = npr.uniform(0, 1)
                if prob_drop < 0.7:
                    lines_normal.append(line)

with gzip.open(os.path.join(path_unified_cancer, "cancer_junctions.gz"), "wt") as f:
    for line in lines_cancer:
        f.writelines(line)

with gzip.open(os.path.join(path_normals, "normal_junctions_drop.gz"), "wt") as f:
    for line in lines_normal:
        f.writelines(line)


lines_expression = []
with gzip.open(os.path.join(path_normals, "normal_junctions_drop.gz"), "rt") as f:
    lines = f.readlines()
    kmer = lines[1:]
    lines_expression.append(lines[0])
    for line in kmer:
        #Set all the expression to 0 with a prob of 0.2
        prob_drop = npr.uniform(0, 1)
        if prob_drop < 0.2:
            line = line.split()
            #Fill the line with zeros
            line[5:] = np.zeros(len(line[5:])).astype(str)
            line = "\t".join(line) + "\n"
            lines_expression.append(line)

        #Now we need to change the expression level of the rest of the kmers to a number between 0 and 150, integers only
        else:
            line = line.split()
            line[5:] = npr.randint(0,105, len(line[5:])).astype(str)
            line = "\t".join(line) + "\n"
            lines_expression.append(line)

with gzip.open(os.path.join(path_normals, "normal_junctions.gz"), "wt") as f:
    for line in lines_expression:
        f.writelines(line)







