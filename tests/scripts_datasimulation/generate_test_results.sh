#!/bin/bash

set -x
set -e

testname="simulated_Ipp"

basedir=$(dirname $(pwd) )

datadir=${basedir}/data_simulated/${testname}

cd ..
### merging splice graphs

outdir=${datadir}/results_merged
mkdir -p $outdir
python -m spladder.spladder build -o ${outdir} -a ${datadir}/${testname}.gtf --bams ${datadir}/align/${testname}_1_sample1.bam,${datadir}/align/${testname}_1_sample2.bam,${datadir}/align/${testname}_1_sample3.bam,${datadir}/align/${testname}_1_sample4.bam,${datadir}/align/${testname}_1_sample5.bam,${datadir}/align/${testname}_1_sample6.bam,${datadir}/align/${testname}_1_sample7.bam,${datadir}/align/${testname}_1_sample8.bam,${datadir}/align/${testname}_1_sample9.bam,${datadir}/align/${testname}_1_sample10.bam,${datadir}/align/${testname}_1_sample11.bam,${datadir}/align/${testname}_1_sample12.bam,${datadir}/align/${testname}_1_sample13.bam,${datadir}/align/${testname}_1_sample14.bam,${datadir}/align/${testname}_1_sample15.bam,${datadir}/align/${testname}_1_sample16.bam,${datadir}/align/${testname}_1_sample17.bam,${datadir}/align/${testname}_1_sample18.bam,${datadir}/align/${testname}_1_sample19.bam,${datadir}/align/${testname}_1_sample20.bam --extract-ase --event-types exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip --readlen 15 --output-conf-icgc --output-txt --output-txt-conf --output-gff3 --output-struc --output-struc-conf --output-bed --output-conf-bed --output-conf-tcga
