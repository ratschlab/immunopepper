# ImmunoPepper
ImmunoPepper is a software tool that takes a splicing graph (possibly derived from an RNA-Seq
samples) as input and generates the set of all theoretically peptide sequences (or kmers) through
direct translation of all walks along the graph. This peptide set can be personalized with germline
and somatic variants and takes all exon-exon junctions present in the splicing graph (even ones not
part of the reference annotation but present in the given RNA-Seq sample) into account. The
comprehensive set of peptides can be used subsequently for further downstream analyses such as
domain annotation or computational immunology.

## Get Started

### Installation

It is recommended to setup a separate virtual or conda environment.
The basic ImmunoPepper package can be installed via `pip`:

```
pip install immunopepper
```

Alternatively, ImmunoPepper can also be installed from source using:
```
pip install -r requirements.txt -r requirements_dev.txt
make install
```

After installation, please consult the help screen for further usage options:
```
immunopepper -h
```
### Prerequisites
ImmunoPepper takes a splicing graph as input. This splicing graph has to be generated using the
SplAdder pipeline. Further information about SplAdder is available on its [GitHub
page](https://github.com/ratschlab/spladder) or the [Online
documentation](https://spladder.readthedocs.io/en/latest/).

## Basic workflow
The software has four basic working modes:

1. `build`: Core part of ImmunoPepper. Traverses the input splice graph and generates all possible
peptides/kmers.
2. `make_bg`: Integrates multiple kmer files (produced via `build`) and generates one background kmer file.
3. `diff`: Takes as input the foreground kmer file and a background kmer file. The output is
contrasting foreground and background, indicating all foreground kmers not present in the
background.
4. `filter`: Apply different filter mechanisms to a given kmer file.

### Mode `build`
The following parameters are *mandatory*:
- `--samples`: input sample names; can specify more than one sample. (Example: 'sample1 sample2')
- `--output-dir`: output directory
- `--ann-path`:annotation file, accepted file formats: *.gtf*, *.gff* and *.gff3*
- `--splice-path`: path of the input SplAdder splice graph
- `--ref-path`: reference genome file in FASTA format

The following parameters are *optional*:
- `--mutation-mode`: mutation mode; choose from {*ref*, *germline*, *somatic* and *somatic_germline*}, default mode *ref*.
- `--kmer`: length of the kmers for kmer ouput. Default value is `0`, which will output full peptides instead of kmers. A recommended kmer length is 9.
- `--disable-concat`: Turns off the generation of kmers from combinations of more than 2 exons (kmers generated from combinations of short exons might be missed)
- `--germline`: germline mutation file path. Mandatory argument if the mutation mode is `germline` or `somatic_and_germline`.
- `--somatic`: somatic mutation file path. Mandatory argument if the mutation mode is `somatic` or `somatic_and_germline`.
- `--use-mut-pickle`: Summarize mutation information in a pickle file and re-use if existing (saves the time processing the original mutation files).
- `--count-path`: path to splice graph count file
- `--gtex-junction-path`: path to intron whitelist. Accept the hdf5 file format. Will accept tsv format in the future.
- `--compressed`: compress the output files using gzip

Example command line (replace `ref` with `germline` to consider mutation information)
```
immunopepper build \
--output-dir tests/test1/current_output_pos \
--ann-path tests/test1/data/test1pos.gtf \
--ref-path tests/test1/data/test1pos.fa \
--splice-path tests/test1/data/posgraph/spladder/genes_graph_conf3.merge_graphs.pickle \
--somatic tests/test1/data/test1pos.maf \
--germline tests/test1/data/test1pos.vcf \
--samples test1pos \
--mutation-mode ref \
--kmer 4 \
--disable-concat \
--count-path  tests/test1/data/posgraph/spladder/genes_graph_conf3.merge_graphs.count.hdf5 \
--gtex-junction-path tests/test1/data/posgraph/spladder/genes_graph_conf3.test1pos.junction.hdf5
```

### Mode `make_bg`
The following parameters are *mandatory*:
- `--kmer-files`: The list of kmer files output by build mode, e.g., 'ref_back_kmer.txt somatic_back_kmer.txt'.
- `--output-file-path`: Output integrated background kmer file path.
- `--output-dir`: Directory to store the log file.

The following parameters are *optional*:
- `--verbose`: Specify the level of output. 0 means zero debug information, 2 means the most detailed information.
- `--compressed`: Compress the output files with gzip.

Example command line:
```
immunopepper make_bg \
--kmer-files tests/test1/current_output_pos/test1pos/ref_back_kmer.txt tests/test1/current_output_pos/test1pos/ref_back_kmer.txt \
--output-dir tests/test1/current_output_pos/ \
--output-file-path tests/test1/current_output_pos/test1pos/uniq_back_kmer.txt \
--verbose 2
```

### Mode `diff`
The following parameters are *mandatory*:
- `--junction-kmer-file`: foreground junction file path generated by `build` mode, e.g., `ref_junction_kmer.txt`
- `--bg-file-path`: background kmer file path. Can be the output of `make_bg` mode or external file. One kmer per line.
- `--output-file-path`: output tsv file path.
- `--output-dir`: directory to store the log file.

The following parameters are *optional*:
- `--verbose`: Specify the verbosity level of output. 0 means zero debug information, 2 means the most detailed information.
- `--compressed`: compress the output files with gzip.

Example command line
```
immunopepper  diff \
--junction-kmer-file tests/test1/current_output_pos/test1pos/ref_junction_kmer.txt \
--bg-file-path tests/test1/current_output_pos/test1pos/uniq_back_kmer.txt \
--verbose 1 \
--output-file-path tests/test1/current_output_pos/test1pos/kmer_result.tsv \
--output-dir tests/test1/current_output_pos \
--remove-bg
```

### Mode `filter`
The following parameters are *mandatory*:
- `--junction-kmer-tsv-file`: The original kmer tsv files. Generated by `build` mode
or by `diff` mode. It should contain field `cross-junction`, `seg-expr` and `junc_expr`.
- `--output-file-path`: Mandatoray argument. Specify the output tsv file path.
- `--output-dir`: Mandatoray argument. Specify the directory to store the log file.

The following parameters are *optional*:
- `--cross-junction`: Only output the cross-junction kmers.
- `--seg-expr`: Only output kmers that have segment expression greater than threshold.
- `--seg-expr-thresh`: Segment expression threshold. Default 0.
- `--junc-expr`: Only output kmers that have junction expression greater than threshold.
- `--junc-expr-thresh`: Junction expression threshold. Default 0.

When providing the `_metadata.tsv.gz` file output by `build` mode, we can have more filters.

- `--meta-file-path`: the meta data file output by `build` mode.
- `--peptide-annotated`: 1 or 0. Filter the kmers based on whether their original kmers appear in background peptide, 0 means keeping the kmers whose original peptide does not show in background peptide. 1 means the opposite.
- `--junction-annotated`: 1 or 0. Filter the kmers based on whether their corresponding junction appear in annotation file, 0 means keeping the kmers whose original junction does not show in annotation file. 1 means the opposite.
- `--has-stop-codon`: 1 or 0. Filter the kmers based on whether their corresponding sequence contains stop codon, 0 means keeping the kmers whose corresponding dna does not contain stop codon. 1 means the opposite.
- `--is-in-junction-list`: 1 or 0. Filter the kmers based on whether their corresponding intron is in the junction whitelist, 0 means keeping the kmers whose corresponding intron id not in the junction whitelist. 1 means the opposite.
- `--is-isolated`: 1 or 0. Filter the kmers based on whether their corresponding peptide comes from single exon or not, 0 means keeping the kmers whose corresponding peptide comes from exon pairs. 1 means the opposite.
- `--verbose`: Specify the level of output. 0 means zero debug information, 2 means the most detailed information.
- `--compressed`: Compress the output files with gzip.

If we have rna-seq data, we can apply more detailed filter. In filter mode, we can also infer the
exact dna positions that can output the kmers.

- `--infer-dna-pos`: if turned on, we will infer the *exact_kmer_pos* for each kmer in the input
junction file `junction-kmer-tsv-file` and write to the `--output-file-path`. Notice the meta data
file `--meta-file-path` is required here.

There will be an additional column in `--output-file-path` compared with the original `junction-kmer-tsv-file`,
*exact_kmer_pos*. It is in the format \[chr\]\_\[strand\]\_\[somatic_var_comb\]\_\[dna_position\]

>*DDAR* ->  *X_+_._11;23*: 4-mer *DDAR* can be translated from chromosome X, positive strand, no somatic mutation
and dna sequence seq\[11:23\]

>*RTHDAR* -> *12_-_130;106_134;124;112;104*: 6-mer *RTHDAR* can be translated from chromosome 12, negative strand, with somatic mutation
in position 130 and 106. and dna sequence complementary(seq\[124:134\]\[::-1\]+seq\[104:112\]\[::-1\]). complementary is to convert
the original base to its complementary base.

Example command line 1 (filter based on metadata):
```
immunopepper filter \
--junction-kmer-tsv-path \
tests/test1/build/pos/test1pos/ref_junction_kmer.txt \
--output-file-path \
tests/test1/current_output_pos/test1pos/kmer_result_filtered.tsv.gz \
--output-dir \
tests/test1/current_output_pos/ \
--meta-file-path \
tests/test1/build/pos/test1pos/ref_metadata.tsv.gz \
--cross-junction \
--junc-expr \
--compressed \
--peptide-annotated 1 \
--junction-annotated 1 \
--has-stop-codon 1 \
--is-isolated 0 \
--is-in-junction-list 1 \
```

Example command line 2 (infer dna position):
```
immunopepper filter \
--junction-kmer-tsv-path tests/test1/current_output_pos/test1pos/ref_junction_kmer.txt \
--output-dir tests/test1/current_output_pos/ \
--output-file-path tests/test1/current_output_pos/test1pos/ref_junction_exact_dna_pos.txt.gz \
--meta-file-path tests/test1/current_output_pos/test1pos/ref_metadata.tsv.gz \
--infer-dna-pos \
--compressed \
```

Once calculate the dna position, you can also use the script `tests/valid_kmer_dna_pos_map.py` to check if the dna position
can generate the exact kmer.

Example command line 3 (validate dna position):
```
cd tests
python tests/valid_kmer_dna_pos_map.py \
--new-junction-file tests/test1/current_output_pos/test1pos/ref_junction_exact_dna_pos.txt.gz \
--sample test1pos \
--mutation-mode ref \
--ref-path tests/test1/data/test1pos.fa \
--somatic tests/test1/data/test1pos.maf \
--germline tests/test1/data/test1pos.vcf \
```
## post-processing guidlines
For further filtering, the user can use the predicted kmers as input for MHC-binding prediction or
use MS databases for further confirmation.

### MHC-Binding
One option for MHC binding prediction is [NetMHC](http://www.cbs.dtu.dk/services/NetMHC/). Using the
predicted kmers as input, `NetMHC` predicts a peptide-MHC class 1 binding score for each sequence
using a neural network.

### Mass spectrometry
Mass spectrometry data can provide further evidence for the presence of a predicted peptide. There
exist several tools for searching a peptide sequence against a MS database, for instance
[OpenMS](https://www.openms.de).

## Output files
There are 5 files for the `build` mode. `mut_mode` refers to `ref`, `somatic`,  `germline` and `somatic_and_germline`.
- **\[mut_mode\]_back_peptides.fa**: Peptides translated from annotation transcripts. Two lines for one output. The first
line is the transcript ID and the second line is the result peptide.
- **\[mut_mode\]_back_kmer.txt**: kmers generated from **\[mut_mode\]_back_peptides.fa**. There are four columns: \[ *kmer*, *gene_name*, *seg_expr*, *is_crossjunction*\]. The first column is
the result kmer, the second column is the transcript ID, the third column is the average segment expression and the final column is the
flag indicating if the kmer is junction kmer. The final column is False for *all* rows in this file.
- **\[mut_mode\]_peptides.fa**: Peptides translated from traversing splicegraph. Two lines for one output. The first
line is the output ID and the second line is the result peptide.
- **\[mut_mode\]_junction_kmer**.txt: kmers generated from **\[mut_mode\]_peptides.fa**. In addition to the same four columns in **\[mut_mode\]_back_kmer.txt**, there is one more
column in this file. *junction_expr*, refers to the junction counts for those kmers that span over
exon junction. For those with *junction_expr* > 0, the flag `is_crossjunction` is True.
- **\[mut_mode\]_metadata.tsv.gz**: Contain details for every junction pairs.

Detail explanation for columns in **\[mut_mode\]_metadata.tsv.gz**
- **output_id**: In the format of \[gene_nama\]:\[first vertex\]_\[second vertex\]:\[somatic variant combination id\]:\[translation start position\]. Like `GENE1:0_2:0:11`.
`GENE1` is the gene name, `0_2` means this junction consists of vertex 0 and vertex 2. `0` means there is no
somatic mutation or it is the first case of all somatic mutation combination cases. `11` means the translation starts from position 11.
- **read_frame**: int (0,1,2). The number of base left to the next junction pair.
- **gene_name**: str. The name of Gene.
- **gene_chr**: str. The Chromosome id where the gene is located.
- **gene_strand**: str ('+', '_'). The strand of gene.
- **mutation_mode**: str ('ref', 'somatic', 'germline', 'somatic_and_germline'). Mutation mode
- **peptide_annotated**: Boolean. Indicate if the junction peptide also appears in the background peptide.
- **junction_peptided**: Boolean. Indicate if the junction also appear in the input annotation file.
- **has_stop_codon**: Boolean. Indicate if there is stop codon in the junction pair.
- **is_in_junction_list**: Boolean. Indicate if the junction pair appear in the given junction whitelist.
- **is_isolated**: Boolean. Indicate if the output peptide is actually translated from a single exon instead of two.
- **variant_comb**: shows the somatic variantion combination used in this line of output. seperate by ';'
    eg. 5;25 means the somatic mutation of position 5 and 25 take effect in this output.
- **variant_seg-expr**: shows the corresponding expression of segments where the corresponding somatic mutation is in.
    eg. 257.0;123.2 means the segment where the somatic mutation in position 5 is in has counts 257.0
- **modified_exons_coor**: Shows exon coordination. Usually we have 4 number start_v1;stop_v1;start_v2;stop_v2. They
    have already absorb reading frame so you can use the coord directly to generate the same output peptide.
- **original_exons_coord**: Shows the original exon coordination.
- **vertex_idx**: shows the vertex id of the given junction. eg 5,6 means this junction pair consists of the fifth and
    sixth vertex.
- **junction_expr**: float. The expression of the junction.
- **segment_expr**: float. The weighted sum of segment expression. We split the junction into segments and compute the segment
    expression with the length-weighted-sum expression.

The `.meta` file is compressed by default in all time. The user can add `--compressed` option
in the input argument to have other files compressed. It is recommended to output in compressed format because
it can save a lot of storage.

The output file for `make_bg` mode is a text file. Each line is a unique kmer.

The output file for `diff` mode is a text file. There is a header line like **\[mut_mode\]_junction_kmer** but
with one more column `is_neo_flag` to indicate if the kmer also exist in the background kmer file. We can also
remove those kmers that exist in the background files with the option `--remove-bg`.

The output file for `filter` mode is a text file also with header line.

## Example use case on experimetal data
Using real DNA-sequencing data from mouse, we will show how to apply ImmunoPepper to generate all
candidate kmers. In this example, we consider two samples: `ENCSR000BZG` and `ERR2130621`. `ERR2130621` comes
from[An RNA-Seq atlas of gene expression in mouse and rat normal tissues](https://www.nature.com/articles/sdata2017185#data-records). It aims to provide a comprehensive RNA-seq dataset
across the same 13 tissues for normal mouse and rat. `ERR2130621` is the RNA-seq result from one of the house mouse.
`ENCSR000BZG` comes from [project ENCODE](https://www.encodeproject.org/report/?type=Experiment&status=released).
It is also a RNA-seq result from central nervous system tissue of house mouse embryo. They are both taken
from normal tissues. In real applications, we should choose the cancer sample as the foreground one and
subtract the normal background result. However, here we just simply show the basic workflow on
real dataset.

We choose `ENCSR000BZG` as the background sample and `ERR2130621` as the foreground sample. They use the same
splicegraph but have different expression values and individual (personalized) mutations. The mutation file
 is created arbitrarily and have no real implications. Our goal is to generate all kmers unique to `ERR2130621`.

### Generate the mouse data
We first show how to generate the data files required to run immunopepper from the files available on public.
- [Mouse Genomic FASTA](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.p6.genome.fa.gz). Get `GRCm38.p6.genome.fa`
- [Mouse Comprehensive gene annotation](https://www.gencodegenes.org/mouse/release_M23.html). Get `gencode.vM23.annotation.gtf`
- [ENCSR000BZG bam file](https://www.encodeproject.org/files/ENCFF713JNY/@@download/ENCFF713JNY.bam) Get `ENCFF713JNY.bam`. We can rename `ENCFF713JNY.bam` as `ENCSR000BZG.bam` to be in consistence with
the given file.
- [ERR2130621 fasta](https://www.ebi.ac.uk/ena/browser/view/ERR2130621). Get `ERR2130621.fastq`

If the binary alignment map file (.bam) is available, we should start from the original RNA-seq fasta data. We
can using alignment toolbox [STAR](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html).
Then we can use [SAMtools](http://samtools.sourceforge.net) to index the alignment. In our case, bam file for sample `ERR2130621` is not available so
we start with the original fasta file.

```
STAR --runThreadN 3 --runMode genomeGenerate --genomeDir genome --genomeFastaFiles  GRCm38.p6.genome.fa --sjdbGTFfile gencode.vM23.annotation.gtf
STAR --runThreadN 3 --genomeDir genome --readFilesIn ERR2130621.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ERR2130621
mv ERR2130621aligned.sortedByCoord.out.bam ERR2130621.bam || TRUE
samtools index  ERR2130621.bam
```
It takes around 6000 seconds to obtain the bam file from `ERR2130621.fastq` with a single thread.

If the binary alignment map file (.bam) is available, we can directly run SplAdder with command to build
 the merge splicegraph.
```
spladder build \
-b ENCSR000BZG.bam ERR2130621.bam\
-o spladder_out \
-a gencode.vM23.annotation.gtf \
-v \
-c 3 \
--merge-strat merge_graphs \
--quantify-graph \
--no-extract-ase \
--no-primary-only \
--no-insert-ir \
--no-insert-es \
--no-insert-ni \
--no-remove-se \
--no-validate-sg \
```
The output of this command line is in `spladder_out/spladder/genes_graph_conf3.ENCFF713JNY.pickle` and the count
file `graph/spladder/genes_graph_conf3.merge_graphs.pickle` and `graph/spladder/genes_graph_conf3.merge_graphs.count.hdf5`.
It takes around 7000 seconds for creating the merged splicegraph using a single thread.

The mutation files *ImmunoPepper_usecase.maf* and *ImmunoPepper_usecase.vcf* are created arbitrarily. We
just choose a position randomly and replace the original base with another one.

### Get the mini-version
Both the `ImmunoPepper_usecase.gtf` and `genome1.fa` are mini-version. The test cases only contain gene `ENSMUSG00000025902.13`
and `ENSMUSG00000025903.14` and they are both in `chr1`. So we can just extract a small part
of the annotation and genome sequence file.

```
grep -E 'ENSMUSG00000025902.13|ENSMUSG00000025903.14' gencode.vM23.annotation.gtf > ImmunoPepper_usecase.gtf
head -n 23252381 GRCm38.p6.genome.fa > genome1.fa
```

The SplAdder output contains more than 50000 genes. In our small usecase, we just contain 2 genes, whose id
is 19 and 33. We can thus extract a mini-version splicegraph and its corresponding expression data, that is *ImmunoPepper_usecase.pickle* and *ImmunoPepper_usecase.count.hdf5*, to accelerate the running process.
The corresponding script is `scripts/create_mini_mouse_data.py`.
You can also use the original data *genes_graph_conf3.merge_graphs.pickle* and *genes_graph_conf3.merge_graphs.count.hdf5* to have a full run.



### Run Immunopepper on mouse data
- Step 1: Use the `build` mode to generate kmers of the two samples in all four mutation modes:
```
cd tests/mouse_usecase
# reference (ref) mode
immunopepper build --mutation-mode ref --samples ENCSR000BZG ERR2130621 --output-dir ImmunoPepper_usecase_out --splice-path ImmunoPepper_usecase.pickle --ann-path ImmunoPepper_usecase.gtf --ref-path genome1.fa --kmer 9 --count-path ImmunoPepper_usecase.count.hdf5
# germline mode
immunopepper build --mutation-mode germline --samples ENCSR000BZG ERR2130621 --output-dir ImmunoPepper_usecase_out --splice-path ImmunoPepper_usecase.pickle --ann-path ImmunoPepper_usecase.gtf --ref-path genome1.fa --kmer 9 --count-path ImmunoPepper_usecase.count.hdf5 --germline ImmunoPepper_usecase.vcf --somatic ImmunoPepper_usecase.maf
# somatic mode
immunopepper build --mutation-mode somatic --samples ENCSR000BZG ERR2130621 --output-dir ImmunoPepper_usecase_out --splice-path ImmunoPepper_usecase.pickle --ann-path ImmunoPepper_usecase.gtf --ref-path genome1.fa --kmer 9 --count-path ImmunoPepper_usecase.count.hdf5 --germline ImmunoPepper_usecase.vcf --somatic ImmunoPepper_usecase.maf
# germline and somatic mode
immunopepper build --mutation-mode somatic_and_germline --samples ENCSR000BZG ERR2130621 --output-dir ImmunoPepper_usecase_out --splice-path ImmunoPepper_usecase.pickle --ann-path ImmunoPepper_usecase.gtf --ref-path genome1.fa --kmer 9 --count-path ImmunoPepper_usecase.count.hdf5 --germline ImmunoPepper_usecase.vcf --somatic ImmunoPepper_usecase.maf
```

- Step 2: Create background kmer set from the output of sample `ENCSR000BZG`.
Since there exist no mutations in sample `ENCSR000BZG`, we only consider its output in reference. In addition,
we only consider kmers that have junction expression larger than 0. We can achieve this using the `filter` mode and
get the file `ref_mode_background_kmer.tsv`. We then use the `make_bg` mode to create the background kmer file. Since
the input is just one file, `make_bg` simply takes the first column and outputs all unique kmers.
```
immunopepper filter --output-dir ImmunoPepper_usecase_out --output-file-path ImmunoPepper_usecase_out/ENCSR000BZG/ref_mode_background_kmer.tsv --junction-kmer-tsv-path ImmunoPepper_usecase_out/ENCSR000BZG/ref_junction_kmer.txt --junc-expr
immunopepper make_bg --kmer-files ImmunoPepper_usecase_out/ENCSR000BZG/ref_mode_background_kmer.tsv --output-dir ImmunoPepper_usecase_out --output-file-path ImmunoPepper_usecase_out/background_kmer.txt
```

- Step 3: Remove the background kmers
After generating the background kmers in Step 2, we can now subtract those kmers from the kmer sets of sample `ERR2130621`. We
can use `diff` for this operation:
```
# contrast ref kmers against background
immunopepper diff --junction-kmer-file  ImmunoPepper_usecase_out/ERR2130621/ref_junction_kmer.txt --bg-file-path ImmunoPepper_usecase_out/background_kmer.txt --output-dir ImmunoPepper_usecase_out --output-file-path ImmunoPepper_usecase_out/ERR2130621/ref_junction_kmer_remove-bg.tsv --remove-bg
# contrast germline kmers against background
immunopepper diff --junction-kmer-file  ImmunoPepper_usecase_out/ERR2130621/germline_junction_kmer.txt --bg-file-path ImmunoPepper_usecase_out/background_kmer.txt --output-dir ImmunoPepper_usecase_out --output-file-path ImmunoPepper_usecase_out/ERR2130621/germline_junction_kmer_remove-bg.tsv --remove-bg
# contrast somatic kmers against background
immunopepper diff --junction-kmer-file  ImmunoPepper_usecase_out/ERR2130621/somatic_junction_kmer.txt --bg-file-path ImmunoPepper_usecase_out/background_kmer.txt --output-dir ImmunoPepper_usecase_out --output-file-path ImmunoPepper_usecase_out/ERR2130621/somatic_junction_kmer_remove-bg.tsv --remove-bg
# contrast somatic/germline kmers against background
immunopepper diff --junction-kmer-file  ImmunoPepper_usecase_out/ERR2130621/somatic_and_germline_junction_kmer.txt --bg-file-path ImmunoPepper_usecase_out/background_kmer.txt --output-dir ImmunoPepper_usecase_out --output-file-path ImmunoPepper_usecase_out/ERR2130621/somatic_and_germline_junction_kmer_remove-bg.tsv --remove-bg
```

- Step 4: Filter
After removing the background kmers in Step 3, we can add more filters to further reduce the number of candidate kmers.
For example, we only consider the kmers that have junction expression larger than 0 as well as a
segment expression value larger than 2. `filter` mode provides filters based on segment expression
and junction expression, based on a user-provided threshold.
```
# filter ref kmers
immunopepper filter --output-dir ImmunoPepper_usecase_out --output-file-path ImmunoPepper_usecase_out/ERR2130621/ref_junction_kmer_remove-bg_filter.tsv --junction-kmer-tsv-path ImmunoPepper_usecase_out/ERR2130621/ref_junction_kmer_remove-bg.tsv --cross-junction --seg-expr --seg-expr-thresh 2
# filter germline kmers
immunopepper filter --output-dir ImmunoPepper_usecase_out --output-file-path ImmunoPepper_usecase_out/ERR2130621/germline_junction_kmer_remove-bg_filter.tsv --junction-kmer-tsv-path ImmunoPepper_usecase_out/ERR2130621/germline_junction_kmer_remove-bg.tsv --cross-junction --seg-expr --seg-expr-thresh 2
# filter somatic kmers
immunopepper filter --output-dir ImmunoPepper_usecase_out --output-file-path ImmunoPepper_usecase_out/ERR2130621/somatic_junction_kmer_remove-bg_filter.tsv --junction-kmer-tsv-path ImmunoPepper_usecase_out/ERR2130621/somatic_junction_kmer_remove-bg.tsv --cross-junction --seg-expr --seg-expr-thresh 2
# filter germline/somatic kmers
immunopepper filter --output-dir ImmunoPepper_usecase_out --output-file-path ImmunoPepper_usecase_out/ERR2130621/somatic_and_germline_junction_kmer_remove-bg_filter.tsv --junction-kmer-tsv-path ImmunoPepper_usecase_out/ERR2130621/somatic_and_germline_junction_kmer_remove-bg.tsv --cross-junction --seg-expr --seg-expr-thresh 2
```

- Step 5: Aggregate
We get the unique kmers of sample `ERR2130621` in four modes. Now we can aggregate all those kmers.
```
tail -n +2 ImmunoPepper_usecase_out/ERR2130621/*_junction_kmer_remove-bg_filter.tsv | cat | grep -v "==>" | cut -f1 | sort |uniq | grep . > neo_kmer.txt
```
## Pratical Tips
- ImmunoPepper requires the sample name are **exactly** the same in the `splice count file` and `mutation file` and the
given option `--samples` should be those samples. Please make necessary changes to the input files so that ImmunoPepper can work as expected.

- `make_bg`, `diff` and `filter` mode accept the output files of ImmunoPepper. However, the user
can also add other external input files.
- `make_bg` assumes the input file has
a header line, separated with *\t* and kmers are in the first column.
- `diff` assumes
the foreground kmer file has a header line and that the background kmer file has the
format as the output file of `make_bg`.
- `filter` assumes the input file has a header
line and with three columns `seg-expr`, `junction_expr` and `is_crossjunction`. It's acceptable if some columns
are missing but the user should not use corresponding filter rules. Otherwise error will happen.

## License
Please see the LICENSE file for more information about license and copyright.
