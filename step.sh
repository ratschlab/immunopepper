#!/bin/bash

set -e

work_dir="/Users/jiayu/PycharmProjects/CBM_RA/projects2018_immunopepper"
spladder_dir="/Users/jiayu/PycharmProjects/CBM_RA/spladder/python"
script_dir=${work_dir}/immunopepper
test_data_dir=${work_dir}/tests


# Step 1and2: make changes on the test1.gtf file and reference sequence
# make some changes on the expression data and Create fq file
python ${script_dir}/create_bam_file.py

# Step 3: Generating genome indexes for STAR
mkdir ${test_data_dir}/data/genome || TRUE
STAR --runThreadN 3 --runMode genomeGenerate --genomeDir ${test_data_dir}/data/genome --genomeFastaFiles ${test_data_dir}/data/test1.fa --sjdbGTFfile ${test_data_dir}/data/test1.gtf --sjdbOverhang 14


## Step 4: Running maping jobs
mkdir ${test_data_dir}/data/align || TRUE
STAR --runThreadN 3 --genomeDir ${test_data_dir}/data/genome --readFilesIn ${test_data_dir}/data/test1_1.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${test_data_dir}/data/align/

## Step 5: Using Samtools to create bai file
samtools index  ${test_data_dir}/data/align/Aligned.sortedByCoord.out.bam

## Step 6: Generating splicegraph and corresponding count file using spladder
rm -rf ${test_data_dir}/data/spladder || TRUE
rm ${test_data_dir}/data/test1.gtf.pickle || TRUE
python ${spladder_dir}/spladder.py  --insert_ir=n --insert_es=n --insert_ni=n --remove_se=n --validate_sg=n -b ${test_data_dir}/data/align/Aligned.sortedByCoord.out.bam -o ${test_data_dir}/data -a ${test_data_dir}/data/test1.gtf -v y -c 3 -M merge_graphs -T n -P n -p n -q y

## Step 7: remove the existing annotation pickle
rm ${test_data_dir}/annotation_preprop.pickle || TRUE

## Step8: run main file to generate result
#python ${script_dir}/main_immuno.py --output_dir tests --ann_path ${test_data_dir}/data/test1.gtf --splice_path ${test_data_dir}/data/spladder/genes_graph_conf3.merge_graphs.pickle --ref_path ${test_data_dir}/data/test1.fa --count_path ${test_data_dir}/data/spladder/genes_graph_conf3.merge_graphs.count.hdf5

