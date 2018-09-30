# immuno
immuno is a script which predicts the pepetide sequency from the RNA sequency splicegraph. You can specify the donor and output using the command line below.


```
python main_immuno.py --donor TCGA-13-1489 --output_dir test --splice_path quick_test_data/sample_gene.pkl --ann_path quick_test_data/small.gtf --ref_path quick_test_data/smallgene34.fa
python exon_filter.py --meta_file test/TCGA-13-1489/ref_meta1.tsv.gz --peptide_file test/TCGA-13-1489/ref_peptide2.fa
```

If you add the expression information, make sure there are `.merge_graphs.pickle` and `.merge_graphs.count.hdf5` files in the same folder.
The github repo provides the splicegraph `sample_gene.pkl` and the annotation information `small.gtf`
but not with the genome sequence `smallgene34.fa` because it is too large.

You should download the genome from other places. We are using the `hg19_hs37d5` version for testing.

Type below command line for more help information.
```
python main_immuno.py -h
```

Please add `modules` folder in [`spladder`](https://github.com/ratschlab/spladder/tree/development/python) repo on the same level with `main_immuno`. The
packages is needed when loading the splicegraph.

### update 30/09/2018
1. Add a new feature to deal with insertion and deletion mutation. To simplify the problem, we
apply the mutation on the list object. eg. To insert `G` behind `A` in position `15`,
we have `dna[15] = `AG`. To delete `G` behind `A` in position `15`, we have
`dna[15:17] = ['A','']`.

2. After considering insertion and deletion, we need to make some changes on the output
filtering part. See `get_exon_dict` in `immuno_filter` model, the element `variant_comb` is
added to the dictionary key.  Since different variantion combination will cause the result to be different
even if they have the same `exon_coord` and `read_frame`.

3. Because of the `read_frame`, the `exon_coord` sometimes will be different from exon
position specified in `splicegraph`. However, the somatic mutation dictionary is attached
to certain exon. In this way, some variant combination is not true. eg. we have variant_comb to be `38;43` while
the `exon_coord` is `40;50;66;74`. It is obvious that only mutation in position `43` takes effect.
This ambiguity will bring trouble to the filtering mechanism. So some changes are made: to compute
the true `variant_comb` in `get_true_variant_comb` function (also in `immuno_filter` module).

4. Build new test files. The new added `test2` case is for the new feature. The original `test1` case
is for the basic requirement.

##### Future work
1. Add k-mer expression data output feature
2. Write good document.


### update 14/09/2018
1. To avoid multimapper case, which might obscure our estimate for the expression count data, the test case
is splitted into positive and negative case. In the following development, all the test case will also
consists two parts.

2. Change the hierarchy of test directory. In future development, we just need create new test case `test2`, `test3`. More pre
cisely, creating `test2pos.fa`, `test2pos.maf`, `test2pos.vcf`, `test2pos.gtf`. Then run `create_bam_file` to generate the corresponding negative case
Follow the next steps stated in `step.sh`. The groundtruth file can be generated automatically in `tests/test2/test2pos` or `tests/test2/test2neg`. Make some
small changes on the `test_end_to_end.py` file. The whole process is done.

>

    .
    ├── ...
    ├── tests                    # root directory for all tests
    │   ├── test1
    │   ├── test2
    │   └── test3
    │       ├── test3pos          # groundtruth for positive case
    │       ├── test3neg          # groundtruth for negative case
    │       └── data              # test data for test3
    │
    │        # Unit tests
    └── ...

3. Move scripts to `./scripts` directory and keep the `./immunopepper` a real packaeg.
4. Calculate all the theoretical value for segment expression in `scripts/validate_count.py`.
5. are incorrect in total 12 segments. Need further check.

##### Future work
Some bash script trick. How to make the `test1`, `test2` flexible variable? like `'insert_a_variable{}'.format('_with_format')`

### update 10/09/2018
1. Add a copy dictionary of `gene.reading_frames` in `immuno_model.py`. Avoid dynamic changes of `gene` object.
2. Add a `constant.py` file and create an important global variable `NOT_EXIST`. It will appear in many places including
**isolated exon case**,  **no mutastion exist case** etc.
3. Decapitalize the namedtuple name.


### update 07/09/2018
1. Split the computation and writing parts of `annotate_gene_opt`. Because of this modification, we can
filter the redundant peptide completely accurate. see `18/06/2018` for more details
2. Use namedtuple to improve the readability of functions.
3. Achieve two new features:
    * calculate the mean segment expression counts. Instead of only replying on the junction expression counts,
    we add one more data to support the expression level of the given exon pair.
    * Generate the libsize level file.

##### Future work
1. Still have some questions about the count file.
    * The count data in the current test case might be inaccurate.
    * why the `strain_id` in `count.h5` file is `Aligned.sortedByCoord.out`, how to change it?
    * How to include more count data in one `count.h5` file?
2. Better to test on the big data to see if something important is ommitted.

### update 24/08/2018
1. Add new test cases. There are 7 transcripts, 7 vertices and 13 output pairs in total in the new test cases. Tests include
    * read_frame propogation
    * somatic and germline mutation
    * filter function
    * isolated exon case
2. When running the test case, fix some bugs:
    * `is_output_redundant`. The current function in `master` branch is wrong because it does not consider the read_frame compare.
    * `get_sub_mut_seq`. Will cause bug when the mutation position is on the boundary. The original one in `master` will cause bug because sometimes the `start_v1` is not
    eqal to the the start position of exon due to read_frame shift.
    * `isolated_peptide_result`. In the negative strand case, the reference should be transformed using `complementary_seq`.
3. Refine the whole process for test case creation in `step.sh`
Put all the steps in the `step.sh` so that it will be quicker to change the test case case in the future.

##### Future work
1. Refactor the code.
    * Split the `immuno_model.py` into two parts for computation and writing. Making the loop of body another function.
    * Using named tuple to avoid complex return value
    * Find a way to avoid making changes on the gene object.
2. Merge into the `fix/envCleanup` branch.
3. Achieve new requirement on the k-mer expression data.


### update 13/08/2018
1. Build the basic test pipline using `pytest` module. It consists of two parts: test for one of the core function `get_sub_mut_dna` and test for all the output file (`peptide.fa`, `metadata.tsv`) in four modes (`ref`, `germline`, `somatic`, `somatic_and_germline`). Type the following command for testing. Remember to include `modules` directory from `spladder` repo in the work directory.
```
python -m pytest test_immunopepper.py
```
2. Generate artifical expression data and corresponding splicegraph count file. The steps are recored in `step.sh` file.

3. Update the test file after including the count file and fix one small bug.


##### Future work
1. More test case for other core functions which may change in the following versions. For example `cross_peptide_result`
2. The artificial `count file` is a little bit different from what was given in the past. Need further confirmation and change if we want to test or use the count file.In the `.count.hdf5`, the `strains` attribute is `Aligned.sortedByCoord.out` instead of `test1`. `gene_names` atrribute shape is (2,1) instead of (2,) so that an additional Index is needed, which is not the case for the previous example file. (See the details In the `parse_gene_metadata_info` function in `immuno_preprocess` module, line 209-211, line 225 and line 231).
3. Need to deal well with `modules` from `spladder` repo. It's a little annoying.

### update 31/07/2018
1. Add the unittest file `test1.fa`, `test1.vcf`, `test1.maf` and the corresponding splicegraph.
Add the function that generates the test file to the `utils` module. That is `create_test_genome`
and `create_full_from_pos`

2. Fix some bugs in different mutation mode after tested on the test file:
    * fix the `sub_dna` creation problem in negative strand case in `get_sub_mut_dna` file.
    * fix the output problem in `germline` mode (which is supposed to output the peptide that not exists in the reference sequence)

##### Future work
1. Create the `.bam` file and based on that create the `.hdf5` encoding the expression information
2. Read relative contents on `pytest`. Combine the current test case with standard python test framework.

### update 21/07/2018
Improve the the speed of mutation mode implementation. In the version before, `mut_seq_dict`
is created, with keys the variant position and values the whole genome sequence. It is a huge waste
because the these different sequences only different in very limited part (just 2 or 3 positions but
we need copy and create new string of about length of over 100 million)

Therefore, instead of copying and pasting the whole genome sequence, we just create a
very small subset of dna string (the range which cover the junction pair). It dramatically
decreases the implementation time.

### update 18/07/2018
1. Add two flags `has_stop_codon` and `is_in_junction_list`. The first new flag indicates if the
output peptide contain stop codon. The underscore `_` symnol is removed in the final output
since the new flag already contains the necessary information. The second new flag is
another way of filtering. From the comparison between the normal sample and cancer sample, we know some
exon pairs exist in gene of both type (normal or cancer) and thus those exon pairs
are less likely to be the key pairs. The new flag is based on the externel file `gtex_junctions.hdf5`

2. Fix a bug caused by different notations. Now all the non-existing value is `.` (
    * the `variant_comb` is `.` if no variant exist in the exon pair;
    * the `seg_expr` and `junction_expr` is `.` if no count info is provided
    * the `junction_annotated` and `is_in_junction_list` is `.` if it is an isolated exon

3. Rearrange the output column. Now all the 5 flags are in the same part.

##### Future work
1. Wait for the suggestion when used in the real application.
2. Need improvement when used in mutation mode. It's both time and memory consuming.

### update 30/06/2018
1. Found some assertion fail.
    * The `cds_stop` does not equal to the `v_stop` (if `strand == -`
    )or `cds_stop` does not equal to the `v_start` (if `strand == +`). Which
    validate further the necessaty of having the reading frame tuple `(cds_start, cds_stop, read_frame)`
    instead of only `read_frame`
    * The isolated cds may not contain stop codon. Which is contrary to our
    basic assumption.

A good example may be the `No.14` (0-based) gene in `splicegraph`. Focus on
the `No.70` vertex. The coordinate is `(50211225,50211389)`. From the corresponding `.gtf` file,
we know there is a cds within that region, with
```
3	HAVANA	CDS	50211318	50211386	.	+	0	gene_id "ENSG00000001617.7"; transcript_id "ENST00000413852.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "SEMA3F"; transcript_type "protein_coding"; transcript_status "NOVEL"; transcript_name "SEMA3F-004"; exon_number 3;  exon_id "ENSE00003521450.1";  level 2; tag "basic"; havana_gene "OTTHUMG00000156806.7"; havana_transcript "OTTHUMT00000345962.1";
```
. It's obvious that it that `50211386` does not equal to `50211389`

Also, the `No.70` vertex has no successive vertice. But it got two read frames.
The above line is one of the read frames.We have
```
seq_dict['3'][50211317:50211386] = ATGTACGTGGGCAGCAAGGACTACGTGCTGTCCCTGGACCTGCACGACATCAACCGCGAGCCCCTCATT
```
The translated peptide is `MYVGSKDYVLSLDLHDINREPLI`, no stop codon.

2. Questions about the new splicegraph. I found the new version splicegraph is different from
the previous one. The most obvious change is the number of vertices: the former splicegraph
has usually around 10 vertices which the new version has more than 400. Quite a lot of
vertices are in overlap. For example, `No.58`, `No.59`, `No.60`, `No.61`, all of those are
connected with `No.70`. The coordinates for the four vertices are
`(5029054,5029254), (50209124, 50209324), (50209426,50209626), (50209485,50209626)`. You can see
they have some overlappings.  `No.67`, `No.68`, `No.69` and `No.70` also very close (share the beginning).
The coordinates are `(50211225,50211259), (_,50211319), (_,50211386), (_,50211389)`

The question is, do we need so many overlapping vertices?

3. Fix the small bug for `exon_filter.py`. Since there are a lot of overlapping vertices
in the new splicegraph. The `exon_filter.py` is more important, which can filter over `40%`
of the peptide. The inner filter module `is_redundant` will not be so useful on this case because most
vertices are in ascending order and as we claim before, in this case, the formal one is redundant but
since it is already written in file so it will not be filtered.

4. Write two logging file `cds_match_log.txt` (`gene_id \t vertex_id` each line)and `isolated_stop` (`gene_name \t vertex_id` each line)

### Update 22/06/18
1. Fix the bug of coordinate mismatch between gene (1-based) and
python (0-based). Now the output peptide on the testset is right.

2. Finalize the filter part on the new structure. Now two ways of
filtering are provided. There is a filter function embedded in `annotate_gene_opt`,
which is `is_output_redundant`. However, since for each output, it
can only be compared with the past output. So there might be some
false negative (the redundant output which should be filtered)
The external way of filtering is the `exon_filter`. Since the global
picture can be seen, the result is 100 percent correct. The reason
why not put the `exon_filter` into `annotate_gene_opt` is that it will
make the code fat (every line of output should be stored). Might
find a better solution in the future.

##### Future work
1. Validate. The previous result can no longer be seen as the fully `ground_truth`
but partly true because there are some missing output. But all the output in the
previous code should also be present in the new output (with a little modification. eg.
stop outputting when meeting with the stop codon). Should find a
better way to validate the current code.

2. Modify the mutation part. The code is right. But each time when
creating new mutated sequence in `apply_somatic_mutation`, a new
sequence string should be created, although just few of them are different,
which means we have huge redundance here. Need to be improved in the future.

### Update 18/06/18
1. Found a small bug in previous code: in the previous code, it is assumed that
the cds always starts at the beginning of vertices (that is `c1=v1` and `c2=v2`),
which actually is not.
2. According to the new requirement, the peptide which contain `stop_codon` should also be outputted until it reaches stop codon.
3. The vertice which has no connection with other vertices but has a cds, called `isolated cds`, should also
be outputted.

Based on the two require, some relative big changes are made on the code.

1. Rewrite the `read_frame` part. In the previous code, the `read_frame` is just the number of bases
that are left to the next vertices. In the new version, the `read_frame` is a tuple
of `(cds_begin_coord, cds_stop_coord,read_frame)`. The first two elements
marks the coding region in the given vertices.
2. Rewrite the vertice pair traverse. First, go through the sorted vertex (ascend if `strand == +` else descend)
If there is a read_frame in that vertex, go to the next inner loop, otherwise read next vertex.
3. Check if the vertex has neighbour, it yes, jump into `cross_peptide_result`.
If there is no stop codon, output the whole peptide result and propagate the `read_frame` to the next vertice.
If there is a stop codon within the vertice pair, output the peptide until it reaches the stop codon. The `read_frame` would
not be propagated. If the stop codon appear at the first vertice (which means no peptide from the second
vertice is outputted), set the flag `is_isolated` to be `True`
4. If the vertex has no neighbor but there is a `cds` in that vertex, it is an `isolated cds`.
The tranlated peptide would also be outputted. `is_isolated = True`, no `read_frame` will be propagated.

##### Future work
1. The result should be validated. When tested on small data, the two tested gene just
output one line of peptide (which means they stop translating at the beginning of coding regions in
reference sequence, which is unbelievable). So far I didn't find the reason.

2. Add the filter part. In the last commit, a filter module was added. Since the
loop was already quite different, a new compatible filter module should be designed.
Can also use `exon_filter_dict`.Or, we can just filter the output based on the
`start_v1/v2, stop_v1/v2`, the four numbers contain all the information needed to output the peptide result.

### Update 12/06/18
1. Validate the exon filter function, add the argument `is_filter` to the interface. When `is_filter` is set to `1`. The exon pair which is the subset of another exon pair will not be processes and the corresponsing peptide will not be output
2.  Build two dictionary for comparing the ground truth output and the current output.

##### Future work
1. Dealing with the `stop codon`. For now the exon pair with `stop codon` will not be translated and outputted, which is not enough for our future analysis.
2. Keep validating.

### Update 28/05/18
1. Add the case where no count file is provided (reference), add new argument `process_num`, which specifies how many genes to be processes in the splicegraph. Default is to process all genes
2. Write bash scripts to generate mini version of vcf and maf file for the future quick validation
3. Keep validating.

### Update 27/05/18
Finish the whole structure.
Below is the flowchar for the core module **annotate\_gene_opt**.

1. Apply mutation information

	* If **.vcf** file is provided, the germline mutation will be encoded. If only **.vcf** is provided (mutation mode: `germline`), the output peptide will be those that are different from the original reference gene.
	* If **.maf** file is provided, the somatic mutation will be encoded. Two dict are created (``som_expr_dict``, ``exon_som_dict``). The first one is used for output the expression info of the segment where the given somatic mutaiton is. The second is used for creating all the variant combination within the given exon pairs.
	* For the `ref` mode, the original peptides are outputted. For the `germline` mode, the peptide which is different from the original peptides are outputted. For the `somatic` mode, the reference sequence is taken as the background sequence and any mutations causing different peptide are outputted. For the `somatic_and_germline` mode. The germline-mutated sequence is taken as the background sequnce.

2. Build two flags for future filtering.
	*	flag0: Indicate if the given peptide output also exist in the expected pepide transcript given by annotation **.gtf** file.
	* flag1: indicate if the junction pair also exist in annotation **.gtf** file
3. Add read phase for each exons
4. loop over each exon pair. Get all variant combnation within the pairs. If any exon combination exists, get the corresponding mutated sequence `mut_seq_dict`.
5. loop over different exon combination and read frames.
6. If the exon pair has no stop codon within the exon pair and the ref_gene is different from mut_gene, or just simplily in `ref` mode. Write the result line on the result file.

The arguments in test.

I try to find the common `donor` in **.h5** and **.maf** file, but failed. So I take the first two donor implied by the **pancan.merged.v0.2.6.PUBLIC.maf** and reduce to a small **.maf** file **mini2.maf** with the mutations that only belong to the two and create a new file **mini_count.h5**which the `strain` are also the two donor (which is actually not). I also only test two genes of over 20,000 genes in the splicegraph.


The complete argument table is shown below.

| Args | input |
|---|---|
|-donor|TCGA-13-1489|
|-output_dir|test|
|-splice_path| ../evaluation\_data\_graph/genes\_graph\_conf3.merge_graphs.validated.pickle |
|-ann_path|/cluster/project/raetsch/lab/03/phrt_immuno/annotation/gencode.v19.annotation.hs37d5\_chr.spladder.gtf|
|-ref_path|../genome/hg19_hs37d5/genome.fa|
|-count_path|../evaluation_data/graph_metadata/genes_graph_conf3.merge_graphs.validated.count.h5|
|-maf|mini.maf|
|-vcf|mini_vcf.h5|

##### Some statistics
1. loading the whole **count.h5** file is time-consuming, around `400` seconds on average. Loading the graph and preprocess the splicegraph are also time-consuming, around `200s` for each.  Loading the whole **.maf** file is even much more time-consuming. So two mini-versioned count and maf file are created.
2. It also need a lot of memory. Usually I 8 CPU(1024 Mb memory for each )to load all these data when submmiting the work via bsub. Four CPUs are not enough.
3. If in `ref` mode, processing each gene costs less than 0.1 seconds.
but for the `somatic` mode, it takes much more time. (from 20s to 80s). Need further analysis to see where the time is spent.

##### Need help
Find real matched **.h5** and **.maf**, **.vcf** file so that I can see how it is real going with the mutation part (currently since no mutation lies within any exon pair, the program never goes into that branch)

##### Future work
1. Tested on more proper dataset
2. Future possible optimization
3. Add more processing message (like some error messages).Currently there is some assumption on the input. Eg. must provide the count file.

### update 05/05/2018
1. Refactor the code so it can be used in four scenarios:`ref`, `germline`, `somatic` and `somatic_and_germline`

    * `ref`: Ouput the peptide on the given splicegraph and reference gene.
    * `germline`: code the germline variant info from **.vcf** file and only
    output the peptide which is different from `ref`
    * `somatic`: code the somatic variants info from **.maf** file and only
    output the peptide which is different from `ref`
    * `somatic_and_germline`: encode the variants info from both **.vcf** file nad **.maf** file.
    using the germline-encoded sequence as the background seq. Only output
    the peptide which is different from what are outputed from background seq.
2. Design a very simple test case: 2 variants in **.maf** and 2
in **.vcf**
3. Create a new immuno module named `immuno_module` which include the core
function `annotate_gene_opt`. It is because we may need do a lot of complicated
loop/if-else control in the future (like dealing with deletion/insertion)
It's better to have a seperate file for the core module.

The code has passed my two self-designed variant file on all four
scenarios.

##### Future work
1. Apply on the real data and improve
2. Might be optimized in terms of speed and memory consumption
3. Deal with more complicated cases (deletion and insertion)

### update 23/04/2018
1. Ouput the junction expression data, which is in the **count.hdf5** and modified by the coefficient in the **lib_size** data
2. Remove some unnecessary argument. Eg. cluster option and load_expression_data.
3. Remove the **get_cfg** module and leave that part directly in the main part.


##### Future work
1. Create new flag based on the expression data
2. Clarify the input data format and keep refactor.

### update 16/04/2018
#### results


1. modularize the code.The details for all modules are shown below:
    * **immuno_preprocess.py'**: include all the preprocess functions (eg: splicegraph,
annotation file, reference gene file).
    * **immuno_filter.py**: include the flags for peptide filtering. Currently two flags.
    * **immuno_mutation.py**: include the function dealing with mutation gene process (to be modified in the future).
    * **inmuno_print.py**: all the print functions might be needed in the project.
    * **utils.py**: all the auxiliary functions needed for output peptides (eg: dna translate, complementary seq)
    * **main_immuno.py**: main functions. Dealing with the input parsing, iterating on the splicegene and call the core function **ann\_gene\_opt**.
    * **get_cfg.py**: parsing the command line input and get the configuration parameters.

2. add more arguments to the main file:

	* **-ann_path**, **splice_path**, **ref_path**. The three arguments can be used to specify the path of the needed input (might have **metagraph_path** in the future). All the path are regarded as the subdirectory of **work_dir**.
	* **-load\_expression_data**. Currently not used in the main function because all the expression annotation data are commented. But it will be used in the future.
	* **quick_test**. Since running on the whole dataset is time-consuming. Quick test mode can help quick validate the code.
	* **cluster_option**. This code will usually run on two different clusters, which result two slightly different settings. Use this option to specify which option **[G/J]**. It will free us from change the code.


##### Future work
1. Start combine the expression data and mutation data in the flags.
2. Keep modularizing.Put all the parameter in the `args`.

### update 08/04/2018
##### results

1. output 3 files.
	* **meta.tsv**: the meta data of output like the gene name, two flags, etc.
	* **peptide.fa**: the output peptide.
	* **no\_output_gene.txt**: record all the genes that are not outputed at all(eg: contain stop codon at very beginning exons)
2. build 2 flags.
	* **flag_background**: indicate whether the output peptide appear in the background peptide.
	* **flag_cds**: indicate whether the exon pair also appear in the cds_list given by *.gtf* annotation file.
3.  create **output_id**: the id looks like **gene_id**.**output_id**. There usually are more than one exon pairs in a gene, so we usually have **0.1**, **4.2**, **108.12**, etc

##### future work
1. consider rare case where the single exon is ignored.
2. modularize the mode.


























