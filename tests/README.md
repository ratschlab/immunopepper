Tutorial folder. This folder contains simulated data and scripts to run the tutorial.

### scripts_datasimulation:

This folder contains the scripts to simulate the data:

- simulate_abundances.py: This script generates the abundances file.

**Input**: The input to this code is simulated_Ipp_blocks.tsv. It is called by using "python simulate_abundances.py {path}/simulated_Ipp_blocks.tsv"
**Output**: simulated_Ipp_abundances.tsv

- generate_testdata_events.sh: This is the main script to generate the data. It calls the script create_event_testcases.py

Within this script, the user should set the prefix for the files (annotation, fastq files...). After calling create_event_testcases.py, the script generates the STAR index for the genome
and aligns the different fastq files against the simulated genome.

**Input**: The script is just called by doing bash generate_testdata_events.sh. However, one should have in the "data_simulated" folder the following files: simulated_Ipp_abundances.tsv and simulated_Ipp_blocks.tsv. The abundances files are generated with the code: "simulate_abundances.py"
**Output**: simulated genome, STAR genome index, annotation file, annotation file for spladder, simulated fastq files and alignment files. The data is stored under the folder "data_simulated/simulated_Ipp"

- create_event_testcases.py: This script generates the genome, the annotation file and the simulated fastq files. The user can choose the length of the reference genome. The annotation file contains trasncript, exon, gene and CDS information.

- generate_test_results.sh: This script calls spladder build mode to generate the count files and the merged splice graph for the input files.
-
#TODO: try to find out how the block file is generated! Or maybe it is just hand made, customized.

#### Steps for data simulation
1. Run simulate_abundances.py to get the abundance file from the block file.
2. Run generate_testdata_events.sh to generate the genome, annotation file and fastq files, as well as the alignment files.
3. Run generate_test_results.sh to call spladder build mode. We will obtain the count files and the merged splice graph for the input files.
4. Create a mutation pileup in bcf format by doing (in command line): samtools mpileup -g -f genome.fa -b alignment_list.txt -o variants.vcf. The list contains a file per line with the path to the alignment files.
5. Convert the bcf to vcf format: bcftools convert -O v -o variants.vcf variants.bcf
6. Run add_mutations.py to generate germline/somatic variations. Run the code changing the path to decide whether we are generating germline or somatic mutations.

The files for these steps are saved in the folder */data_simulated/simulated_Ipp*.



[//]: # (The code simulate_abundances.py takes as an input the file with the blocks.tsv. This file looks like:)

[//]: # (gene1   t1      1       +       51:200,301:400,601:750)

[//]: # (gene1   t2      1       +       51:200,601:750)

[//]: # ()
[//]: # (In the code we choose the number of samples &#40;N=20&#41;, and we simulate the libray sizes: [2.  3.5 4.5 3.5 4.5 4.  2.  3.5 1.  1.5 3.  3.  0.5 3.  0.5 2.  1.  4., 2.5 1. ])

[//]: # (Then, we simulate the total gene expression abundance &#40;total transcript counts per gene&#41;, which is a random number between 20 and 50.)

[//]: # (Finally, the abundance for each gene in each sample is generated. In our example I believe the gene1 is chosen to be differentially expressed)

[//]: # (The output is a file called abundances.tsv, which looks like:)

[//]: # (gene1   t1      48      37      69      49      103     62      44      71      10      27      7       7       2       12      1       8       4       21      12      4)

[//]: # (gene1   t2      5       10      19      13      18      24      11      14      5       8       36      29      8       69      5       24      23      42      59      12)

[//]: # ()
[//]: # (create_event_testcases.py: This code creates the simulated genome of lenght L. This length can be customized by the user. It also creates the gtf file. This gtf file is already modified so that it includes transcript name and CDS information. The cds is created for each transcript, by taking 10 positions after the start to 20 positions before the end. Done in the function &#40;gen_annotation&#41;)

[//]: # (It also creates fq files for each of the samples, generating reads taken from the genome, with a read number that depends of the abundances selected in _abundances.tsv)

[//]: # ()
[//]: # (In order to run this code, one can use generate_testdata_events.sh. It will automatically run the code. Then, it also generates genome indexes and aligns the fq files against the reference simulated genome, using star.)

[//]: # ()
[//]: # ()
[//]: # ()
