import pandas as pd
import gzip
import glob
#First we check the three normal conditions:
#1. No kmers with any expression in more than 15 samples are selected
#2. No kmers with expression bigger than 100 in one sample are selected
#3. No kmers with no expression in any sample selected

path_data ='../../../immunopepper_usecase/filter_case/ref_graph_kmer_NormalExpr/normal_junctions.gz'
path_results_s = '../../../immunopepper_usecase/filter_case/interm_normals_combiExprCohortLim0.0Across1.tsv.gz/part-00000.gz'
path_results_e = '../../../immunopepper_usecase/filter_case/interm_normals_combiExprCohortLim100.0Across1.tsv.gz/part-00000.gz'
path_results_annot = '../../../immunopepper_usecase/filter_case/kmers_derived_solely_from_annotation.tsv.gz/part-00000-2e497c34-7e0e-4a37-88f8-8610b90b771b-c000.csv.gz'
final_normal_results = '../../../immunopepper_usecase/filter_case/final_normal/*part*'
#Read the data to dataframes

with gzip.open(path_data, 'rt') as f:
    df_data = pd.read_csv(f, sep='\t')

with gzip.open(path_results_s, 'rt') as f:
    df_results_s = pd.read_csv(f, sep='\t', header=None)

with gzip.open(path_results_e, 'rt') as f:
    df_results_e = pd.read_csv(f, sep='\t', header=None)

with gzip.open(path_results_annot, 'rt') as f:
    df_results_annot = pd.read_csv(f, sep='\t', header=0)

normal_results= glob.glob(final_normal_results)
normal_res = pd.DataFrame()
for res in normal_results:
    with gzip.open(res, 'rt') as f:
        kmer_list = (pd.read_csv(f, sep='\t', header=0))
        normal_res = pd.concat([normal_res,kmer_list])

#Check condition 1. In this case we need to filter the data to see which kmers have expression in more than 15 samples

df_results_s_filt = df_results_s[df_results_s.iloc[:, 1]>= 15]

# Now I need to filter the data to see which kmers have expression in more than 15 samples.

condition1 = df_data.iloc[:, 5:] > 0.0
condition1 = condition1.sum(axis=1) > 15

# Concatenate my condition1 with the kmers
cond1_final = pd.concat([df_data["kmer"], condition1], axis=1)
kmers_cond1_normals = cond1_final[cond1_final.iloc[:, 1] == True]["kmer"].reset_index(drop=True)

cond1= (df_results_s_filt.iloc[:,0].rename("kmer") == kmers_cond1_normals).all()

if cond1:
   print('Normal samples condition 1 is met. The fitering on sample number is correct')

#Check condition 2. If there is any kmer with expression bigger than 100 we select it as normal

kmers_result2 = df_results_e.iloc[:,0].rename("kmer")

condition2 = df_data.iloc[:, 5:] >= 100.0
condition2 = condition2.any(axis=1)

cond2_final = pd.concat([df_data["kmer"], condition2], axis=1)
kmers_cond2_normals = cond2_final[cond2_final.iloc[:, 1] == True]["kmer"].reset_index(drop=True)

cond2= (kmers_result2 == kmers_cond2_normals).all()

if cond2:
    print('Normal samples condition 2 is met. The fitering on expression is correct')

#Check condition 3. If the annot kmers were correctly identified and removed

df_data_annot = df_data.iloc[:, 7:] == 0.0
df_data_annot = df_data_annot.all(axis=1)
df_kmers_annot = df_data[df_data_annot]["kmer"].reset_index(drop=True)

cond3_identified = (df_results_annot["kmer"] == df_kmers_annot).all()

#Now check if they are in the results. They are in the intermediate files

if (df_kmers_annot.isin(normal_res).any()==False) and cond3_identified:
    print('Normal samples condition 3 is met. Annotation only kmers are identified and removed')
    cond3 = True


if cond1 and cond2 and cond3:
    print('The filtering for the normals is correct')

#Now we check the cancer conditions:
#1. No kmers with junctAnnotated True are selected
#2. Sample filter: Check if kmers with expression >= 20 in sample 3 are selected
#3. Cohort filter: Check if kmers with expression >= 110 in 2 samples at least are selected (excluding sample 3). I can take the simulated file 110 across 1 and pass the threshold of across 2.
#With this I have the kmers resulting of step 3
#Do the union of the kmers resulting of step 2 and step 3, do differential filtering and that is the output of the program. I will save an intermediate file with the kmers resulting of step 2 for the purpose of this check.

#Load data
path_data_cancer = '../../../immunopepper_usecase/filter_case/ref_graph_kmer_CancerExpr/cancer_junctions.gz'
path_results_sample_cancer= '../../../immunopepper_usecase/filter_case/condition2/part-00000-f4e75feb-4886-4f2a-8315-f4aa35d4907c-c000.csv.gz'
path_results_cohort_cancer = '../../../immunopepper_usecase/filter_case/interm_cancer_ref_combiExprCohortLim110.0Across1ExceptsimulatedIpp1sample3.tsv.gz/part-00000.gz'

with gzip.open(path_data_cancer, 'rt') as f:
    df_data_cancer = pd.read_csv(f, sep='\t')

with gzip.open(path_results_sample_cancer, 'rt') as f:
    df_results_sample = pd.read_csv(f, sep='\t', header=0)

with gzip.open(path_results_cohort_cancer, 'rt') as f:
    df_results_cohort = pd.read_csv(f, sep='\t', header=None)

#Check condition 1

df_neo = df_data_cancer[df_data_cancer['junctionAnnotated'] == False]
kmers_neo = df_neo["kmer"]

cond1_a = (kmers_neo == df_results_sample.iloc[0,:].reindex(kmers_neo.index).rename("kmer")).any()
cond1_b = (kmers_neo == df_results_cohort.iloc[0,:].reindex(kmers_neo.index).rename("kmer")).any()

if not (cond1_b and cond1_a):
    print('Cancer condition 1 is met. The filtering on JunctionAnnotated is applied properly')

#Check condition 2

df_results2 = df_results_sample["kmer"].reset_index(drop=True)
condition2 = df_neo.iloc[:, 7] >= 20.0
condition2_kmers = df_neo.loc[condition2, "kmer"].reset_index(drop=True)
cond2 = (df_results2 == condition2_kmers).all()
#
if cond2:
    print('Cancer condition 2 is met. The filtering on sample threshold is correct.')

#Check condition 3

df_results_cohort_filt = df_results_cohort[df_results_cohort.iloc[:, 1]>= 2]
df_results_cohort_kmers = df_results_cohort_filt.iloc[:,0].rename("kmer")
#Now I need to filter the data to see which kmers have expression >=110 in more than 2 samples.

condition3 = df_neo.iloc[:, 5:] >= 110.0
condition3 = condition3.sum(axis=1) >= 2

# Concatenate my condition3 with the kmers
cond3_final = pd.concat([df_neo["kmer"], condition3], axis=1)
kmers_cond3 = cond3_final[cond3_final.iloc[:, 1] == True]["kmer"].reset_index(drop=True)

cond3= (df_results_cohort_filt.iloc[:,0].rename("kmer") == kmers_cond3).all()

if cond3:
    print('Cancer condition 3 is met. The fitering on cohort is correct')

#Check for differential filtering

# #For cancer kmers I will do the intersection between kmers_cond3 and condition2_kmers
# df_cond2 = pd.DataFrame(condition2_kmers)
# df_cond3 = pd.DataFrame(kmers_cond3)
# df_cancer = pd.merge(df_cond2, df_cond3, how='inner', on=['kmer'])
#
# #TODO: how to do the script for normals? The filtering with the output normals and cancer works. I need to find out why the combination of normals does not work
#
# #For the normals I need to join the kmers resulting of condition 1 and condition 2
# df_normal_cond1 = pd.DataFrame(kmers_cond1_normals)
# df_normal_cond2 = pd.DataFrame(kmers_cond2_normals)
#
# df_normal = pd.merge(df_normal_cond1, df_normal_cond2, how='outer', on=['kmer']).reset_index(drop=True)
# #Remove the kmers found only in the annotation
# df_normal = df_normal[~df_normal.kmer.isin(df_kmers_annot)]
#
# df_common = pd.merge(df_normal, df_cancer, how='inner', on=['kmer'])
# #Get the entries of df_cancer that are not present in df_common
# df_neopeptides = df_cancer[~df_cancer.kmer.isin(df_common.kmer)]
# #print(df_neopeptides)
#
# df_neo_res = df_cancer[~df_cancer.kmer.isin(normal_res)]
# print(df_neo_res)
