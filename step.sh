#!/bin/bash

set -e

work_dir="/Users/jiayu/PycharmProjects/CBM_RA/projects2018_immunopepper"
spladder_dir="/Users/jiayu/PycharmProjects/CBM_RA/spladder/python"
testname='test1'

script_dir=${work_dir}/scripts
test_data_dir=${work_dir}/tests/$testname


# Step 1and2: make changes on the test1.gtf file and reference sequence
# make some changes on the expression data and Create fq file
#python ${script_dir}/create_bam_file.py

# Step 3: Generating genome indexes for STAR
#mkdir -p ${test_data_dir}/data/genome/pos || TRUE
#mkdir -p ${test_data_dir}/data/genome/neg || TRUE
#STAR --runThreadN 3 --runMode genomeGenerate --genomeDir ${test_data_dir}/data/genome/pos --genomeFastaFiles ${test_data_dir}/data/test1pos.fa --sjdbGTFfile ${test_data_dir}/data/test1pos.gtf --sjdbOverhang 14 --sjdbInsertSave Basic
#STAR --runThreadN 3 --runMode genomeGenerate --genomeDir ${test_data_dir}/data/genome/neg --genomeFastaFiles ${test_data_dir}/data/test1neg.fa --sjdbGTFfile ${test_data_dir}/data/test1neg.gtf --sjdbOverhang 14 --sjdbInsertSave Basic


## Step 4: Running maping jobs
#mkdir ${test_data_dir}/data/align || TRUE
#STAR --runThreadN 3 --genomeDir ${test_data_dir}/data/genome/pos --readFilesIn ${test_data_dir}/data/test1_1.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${test_data_dir}/data/align/pos
#STAR --runThreadN 3 --genomeDir ${test_data_dir}/data/genome/neg --readFilesIn ${test_data_dir}/data/test1_1.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${test_data_dir}/data/align/neg

# Step 5: Rename the bam file so that we can control the strain id
#mv ${test_data_dir}/data/align/posAligned.sortedByCoord.out.bam ${test_data_dir}/data/align/test1pos.bam || TRUE
#mv ${test_data_dir}/data/align/negAligned.sortedByCoord.out.bam ${test_data_dir}/data/align/test1neg.bam || TRUE

## Step 6: Using Samtools to create bai file
#samtools index  ${test_data_dir}/data/align/test1pos.bam
#samtools index  ${test_data_dir}/data/align/test1neg.bam

## Step 7: Generating splicegraph and corresponding count file using spladder
#rm ${test_data_dir}/data/test1pos.gtf.pickle || TRUE
#rm ${test_data_dir}/data/test1neg.gtf.pickle || TRUE
#rm -rf ${test_data_dir}/data/posgraph || TRUE
#rm -rf ${test_data_dir}/data/neggraph || TRUE
#python ${spladder_dir}/spladder.py  --insert_ir=n --insert_es=n --insert_ni=n --remove_se=n --validate_sg=n -b ${test_data_dir}/data/align/test1pos.bam -o ${test_data_dir}/data/posgraph -a ${test_data_dir}/data/test1pos.gtf -v y -c 3 -M merge_graphs -T n -P n -p n -q y
#python ${spladder_dir}/spladder.py  --insert_ir=n --insert_es=n --insert_ni=n --remove_se=n --validate_sg=n -b ${test_data_dir}/data/align/test1neg.bam -o ${test_data_dir}/data/neggraph -a ${test_data_dir}/data/test1neg.gtf -v y -c 3 -M merge_graphs -T n -P n -p n -q y

## Step8: run main file to generate ground result
# positive case
python ${script_dir}/main_immuno.py --mutation_mode ref --output_dir ${test_data_dir} --ann_path ${test_data_dir}/data/test1pos.gtf --splice_path ${test_data_dir}/data/posgraph/spladder/genes_graph_conf3.merge_graphs.pickle --ref_path ${test_data_dir}/data/test1pos.fa --count_path ${test_data_dir}/data/posgraph/spladder/genes_graph_conf3.merge_graphs.count.hdf5 --samples test1pos --debug --vcf_path ${test_data_dir}/data/test1pos.vcf --maf_path ${test_data_dir}/data/test1pos.maf
python ${script_dir}/main_immuno.py --mutation_mode somatic --output_dir ${test_data_dir} --ann_path ${test_data_dir}/data/test1pos.gtf --splice_path ${test_data_dir}/data/posgraph/spladder/genes_graph_conf3.merge_graphs.pickle --ref_path ${test_data_dir}/data/test1pos.fa --count_path ${test_data_dir}/data/posgraph/spladder/genes_graph_conf3.merge_graphs.count.hdf5 --samples test1pos --debug --vcf_path ${test_data_dir}/data/test1pos.vcf --maf_path ${test_data_dir}/data/test1pos.maf
python ${script_dir}/main_immuno.py --mutation_mode germline --output_dir ${test_data_dir} --ann_path ${test_data_dir}/data/test1pos.gtf --splice_path ${test_data_dir}/data/posgraph/spladder/genes_graph_conf3.merge_graphs.pickle --ref_path ${test_data_dir}/data/test1pos.fa --count_path ${test_data_dir}/data/posgraph/spladder/genes_graph_conf3.merge_graphs.count.hdf5 --samples test1pos --debug --vcf_path ${test_data_dir}/data/test1pos.vcf --maf_path ${test_data_dir}/data/test1pos.maf
python ${script_dir}/main_immuno.py --mutation_mode somatic_and_germline --output_dir ${test_data_dir} --ann_path ${test_data_dir}/data/test1pos.gtf --splice_path ${test_data_dir}/data/posgraph/spladder/genes_graph_conf3.merge_graphs.pickle --ref_path ${test_data_dir}/data/test1pos.fa --count_path ${test_data_dir}/data/posgraph/spladder/genes_graph_conf3.merge_graphs.count.hdf5 --samples test1pos --debug --vcf_path ${test_data_dir}/data/test1pos.vcf --maf_path ${test_data_dir}/data/test1pos.maf

# negative case
python ${script_dir}/main_immuno.py --mutation_mode ref --output_dir ${test_data_dir} --ann_path ${test_data_dir}/data/test1neg.gtf --splice_path ${test_data_dir}/data/neggraph/spladder/genes_graph_conf3.merge_graphs.pickle --ref_path ${test_data_dir}/data/test1neg.fa --count_path ${test_data_dir}/data/neggraph/spladder/genes_graph_conf3.merge_graphs.count.hdf5 --samples test1neg --debug --vcf_path ${test_data_dir}/data/test1neg.vcf --maf_path ${test_data_dir}/data/test1neg.maf
python ${script_dir}/main_immuno.py --mutation_mode somatic --output_dir ${test_data_dir} --ann_path ${test_data_dir}/data/test1neg.gtf --splice_path ${test_data_dir}/data/neggraph/spladder/genes_graph_conf3.merge_graphs.pickle --ref_path ${test_data_dir}/data/test1neg.fa --count_path ${test_data_dir}/data/neggraph/spladder/genes_graph_conf3.merge_graphs.count.hdf5 --samples test1neg --debug --vcf_path ${test_data_dir}/data/test1neg.vcf --maf_path ${test_data_dir}/data/test1neg.maf
python ${script_dir}/main_immuno.py --mutation_mode germline --output_dir ${test_data_dir} --ann_path ${test_data_dir}/data/test1neg.gtf --splice_path ${test_data_dir}/data/neggraph/spladder/genes_graph_conf3.merge_graphs.pickle --ref_path ${test_data_dir}/data/test1neg.fa --count_path ${test_data_dir}/data/neggraph/spladder/genes_graph_conf3.merge_graphs.count.hdf5 --samples test1neg --debug --vcf_path ${test_data_dir}/data/test1neg.vcf --maf_path ${test_data_dir}/data/test1neg.maf
python ${script_dir}/main_immuno.py --mutation_mode somatic_and_germline --output_dir ${test_data_dir} --ann_path ${test_data_dir}/data/test1neg.gtf --splice_path ${test_data_dir}/data/neggraph/spladder/genes_graph_conf3.merge_graphs.pickle --ref_path ${test_data_dir}/data/test1neg.fa --count_path ${test_data_dir}/data/neggraph/spladder/genes_graph_conf3.merge_graphs.count.hdf5 --samples test1neg --debug --vcf_path ${test_data_dir}/data/test1neg.vcf --maf_path ${test_data_dir}/data/test1neg.maf

# step9: unzip gz file and rename to _gt file
python ${script_dir}/generate_gt_file.py
