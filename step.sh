#! /usr/bin/bash

# Step 1: Generating genome indexes for STAR
#STAR --runThreadN 3 --runMode genomeGenerate --genomeDir quick_test_data/genome --genomeFastaFiles quick_test_data/test1.fa --sjdbGTFfile quick_test_data/test1.gtf --sjdbOverhang 14

# Step 2: Running maping jobs
STAR --runThreadN 3 --genomeDir quick_test_data/genome --readFilesIn quick_test_data/test1_3.fq --outSAMtype BAM SortedByCoordinate

# Step 3: Using Samtools to create bai file
samtools index  Aligned.sortedByCoord.out.bam

# Step 4: Generating splicegraph and corresponding count file using spladder
python ../spladder/python/spladder.py  --insert_ir=n --insert_es=n --insert_ni=n --remove_se=n --validate_sg=n -b Aligned.sortedByCoord.out.bam -o quick_test_data -a quick_test_data/test1.gtf -v y -c 3 -M merge_graphs -T n -P n -p n -q y

