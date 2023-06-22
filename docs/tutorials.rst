Tutorials
==========

In this page, a list of tutorials is provided to learn how to use the software. First, the input data used for this tutorial will be described. Then, a tutorial for each mode will be presented.

1. :ref:`input_data`: Description of the input data.
2. :ref:`tutorial_build`: Tutorial for the *build* mode.
3. :ref:`tutorial_samplespecif`: Tutorial for the *samplespecif* mode.
4. :ref:`tutorial_cancerspecif`: Tutorial for the *cancerspecif* mode.
5. :ref:`tutorial_mhcbind`: Tutorial for the *mhcbind* mode.

.. _input_data:

Input data
--------------

The input data used for this tutorial can be found in the folder *'tests/data_simulated/data'*. The data is composed of the following files:

1. **genes_graph_conf3.merge_graphs.pickle:** Splice graphs obtained from SplAdder. The file will contain one splice graph per gene.

    - This file contains information for 9 genes.
    - The object obtained from SplAdder contains data and metadata. However, in this software, only the data is used, which contains different fields. The whole list of the fields can be explored in the SplAdder documentation #TODO: add link. Here, only the relevant fields for this tutorial are shown.

        - **chr:** Chromosome information. In this example all the genes belong to chromosome 1. The chromosomes have the notation "x", x being a valid chromosome number. *Note*: The chromosomes in the reference genome and in the reference annotation file must have the same notation.
        - **name:** Gene ID.
        - **start:** Start coordinate of the genes.
        - **stop:** Stop coordinate of the genes.
        - **strand:** Strand where the genes are placed. It can be either "+" or "-".

    - Therefore, the data contained in the splice graph file for this particular example is the following:


        +------------------+------+-------+------+----------+
        |      chr         | name | start | stop | strand   |
        +==================+======+=======+======+==========+
        |       1          | gene1| 50    | 750  | "+"      |
        +------------------+------+-------+------+----------+
        |       1          | gene3| 1050  | 1600 | "+"      |
        +------------------+------+-------+------+----------+
        |       1          | gene4| 2100  | 2850 | "+"      |
        +------------------+------+-------+------+----------+
        |       1          | gene6| 4050  | 4500 | "+"      |
        +------------------+------+-------+------+----------+
        |       1          | gene7| 4550  | 5000 | "-"      |
        +------------------+------+-------+------+----------+
        |       1          | gene8| 5050  | 5500 | "+"      |
        +-----+------------+------+-------+------+----------+
        |       1          | gene9| 5550  | 6000 | "-"      |
        +------------------+------+-------+------+----------+
        |       1          | gene10| 6050 | 7000 | "+"      |
        +------------------+------+-------+------+----------+
        |       1          | gene15| 9050 | 9450 | "+"      |
        +------------------+------+-------+------+----------+

2. **genes_graph_conf3.merge_graphs.count.hdf5:** Count file obtained from SplAdder.

    - It contains a field *"strains"* that stores the sample names of the processed data. In this case the samples are: "simulated_Ipp_1_samplex", with x being a number from 1 to 20.
    - These sample names are the ones accepted by the software wherever a sample id is required. For example, if one wants to set a `--mutation sample` in the command line, the sample name must be one of the "strains" names.

3. **genome.fa:** Reference genome. It is a fasta file containing the reference genome. For this example, the genome is generated randomly, so it does not corresponds to a biologically accurate genome. The chromosome names must be the same as in the annotation file and in the splice graph file.
4. **genome.fa.fai:** Genome index. The file is useful to access different positions of the genome faster.
5. **simulated_Ipp.gtf:** Reference annotation file. It is a gtf file containing the reference annotation used for this example. The annotation file needs to have information about the gene, gene transcripts, gene exons and gene CDS.
6. **variants_somatic.vcf:** VCF file containing the somatic variants. These somatic variants are randomly generated for the purpose of the example.
7. **variants_germline.vcf:** VCF file containing the germline variants. These germline variants are randomly generated for the purpose of the example.

.. _tutorial_build:

Running build mode
----------------------

Build mode *without* mutations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The main mode of immunopepper is the *build* mode. This mode traverses the splice graph and generates all possible peptides/kmers. It is the base of the other modes. The build mode accepts many different arguments, which add different functionalities to the mode.

First of all, we will show an example with the basic commands, without any additional functionalities such as mutations. The command used for this part of the tutorial is:

.. code-block::

    immunopepper  build --output-dir immunopepper_usecase/ --ann-path  immunopepper/tests/data_simulated/data/simulated_Ipp.gtf --splice-path  immunopepper/tests/data_simulated/data/genes_graph_conf3.merge_graphs.pickle --ref-path  immunopepper/tests/data_simulated/data/genome.fa --kmer 9 --count-path immunopepper/tests/data_simulated/data/genes_graph_conf3.merge_graphs.count.hdf5 --parallel 1 --batch-size 1  --start-id 0 --process-num 0 --output-fasta --verbose 2

By calling this command, the software will generate the possible kmers/peptides for each of the 9 genes in the splice graph. It will take into account the reference genome and the annotation file, and it will generate an output for both the background and foreground peptides. The command is run on the :ref:`input_data` described in the section above. Moreover, the output directory is set to a folder called *immunopepper_usecase*, located on the directory where the command is executed. The kmer length is set to 9, as it is a common kmer length selected in clinical applications.

**Terminal output:**

The output displayed in the command line is the following:

.. code-block:: console

    2023-06-22 12:48:54,100 INFO     Command lineNamespace(output_dir='immunopepper_usecase/', ann_path='immunopepper/tests/data_simulated/data/simulated_Ipp.gtf', splice_path='immunopepper/tests/data_simulated/data/genes_graph_conf3.merge_graphs.pickle', ref_path='immunopepper/tests/data_simulated/data/genome.fa', kmer=9, libsize_extract=False, all_read_frames=False, count_path='immunopepper/tests/data_simulated/data/genes_graph_conf3.merge_graphs.count.hdf5', output_samples=[], heter_code=0, compressed=True, parallel=1, batch_size=1, pickle_samples=[], process_chr=None, complexity_cap=None, genes_interest=None, start_id=0, process_num=0, skip_annotation=False, libsize_path=None, output_fasta=True, force_ref_peptides=False, filter_redundant=False, kmer_database=None, gtex_junction_path=None, disable_concat=False, disable_process_libsize=False, mutation_sample=None, germline='', somatic='', sample_name_map=None, use_mut_pickle=False, verbose=2)
    2023-06-22 12:48:54,100 INFO     >>>>>>>>> Build: Start Preprocessing
    2023-06-22 12:48:54,100 INFO     Building lookup structure ...
    2023-06-22 12:48:54,101 INFO            Time spent: 0.000 seconds
    2023-06-22 12:48:54,102 INFO            Memory usage: 0.159 GB
    2023-06-22 12:48:54,102 INFO     Loading count data ...
    2023-06-22 12:48:54,104 INFO            Time spent: 0.002 seconds
    2023-06-22 12:48:54,104 INFO            Memory usage: 0.160 GB
    2023-06-22 12:48:54,104 INFO     Loading splice graph ...
    2023-06-22 12:48:54,105 INFO            Time spent: 0.000 seconds
    2023-06-22 12:48:54,105 INFO            Memory usage: 0.161 GB
    2023-06-22 12:48:54,105 INFO     Add reading frame to splicegraph ...
    2023-06-22 12:48:54,107 INFO            Time spent: 0.002 seconds
    2023-06-22 12:48:54,107 INFO            Memory usage: 0.161 GB
    2023-06-22 12:48:54,107 INFO     >>>>>>>>> Finish Preprocessing
    2023-06-22 12:48:54,107 INFO     >>>>>>>>> Start traversing splicegraph
    2023-06-22 12:48:54,107 INFO     >>>> Processing output_sample cohort, there are 9 graphs in total
    2023-06-22 12:48:54,108 INFO     Saving results to immunopepper_usecase/cohort_mutNone
    2023-06-22 12:48:54,108 INFO     Not Parallel
    2023-06-22 12:48:54,108 INFO     >>>>>>>>> Start Background processing
    2023-06-22 12:48:54,111 INFO     Saved ref_annot_peptides.fa.gz with 40 lines in 0.0003s
    2023-06-22 12:48:54,111 INFO     Saved ref_annot_kmer.gz with 294 lines in 0.0002s
    2023-06-22 12:48:54,113 DEBUG    ....cohort: annotation graph from batch all/9 processed, max time cost: 0.0, memory cost: 0.16 GB
    2023-06-22 12:48:54,113 INFO     >>>>>>>>> Start Foreground processing
    2023-06-22 12:48:54,175 INFO     Saved gene_expression_detail.gz with 9 lines in 0.0006s
    2023-06-22 12:48:54,176 INFO     Saved ref_sample_peptides.fa.gz with 88 lines in 0.0004s
    2023-06-22 12:48:54,177 INFO     Saved ref_sample_peptides_meta.gz with 44 lines in 0.0005s
    2023-06-22 12:48:54,177 DEBUG    ....cohort: output_sample graph from batch all/9 processed, max time cost: 0.02, memory cost: 0.16 GB
    2023-06-22 12:48:54,188 INFO     Saved library size results to immunopepper_usecase/expression_counts.libsize.tsv

**Output files:**

The output files are saved in the directory *immunopepper_usecase/cohort_mutNone*. The output files are:

1. **ref_annot_peptides.fa.gz:**

    This is a fasta file containing the background peptides for each gene transcript. The file contains a header, that is the transcript id, and the sequence of the corresponding peptide. The name also shows the mutation mode, which in this case is reference. *Note:* As this genome is simulated, there is a higher frequency of stop codons than in nature, that explains the existence of some short peptides.

    This file shows the full transcripts found in the organism as described by the annotation file provided under --ann-path. The sequence of exons for each given transcript is obtained from the annotation file, and the regions corresponding to this exons are taken from the reference genome file provided under `--ref-path`, and translated to create the set of *background* peptides or kmers.

    *Output example:*

    .. code-block::

        fasta
        >gene8.t1
        TSSRTMETLVP
        >gene1.t1
        LQHNSTRSIFWH
        >gene10.t1
        LSLVHPGTRRITKRRRQYPYVIASCQREAGCRGIICS

2. **ref_annot_kmer.gz**: This file contains the kmers of length 9 obtained from the background peptides. The kmers are obtained by passing a sliding window through the peptides contained in ref_annot_peptides.fa.gz. *Note*: For peptides that are shorter than 9 aminoacids, the kmers are not obtained.

    In the example below, one can see the kmers obtained from *gene8.t1*.

    *Output example:*

    .. code-block::

        kmer
        TSSRTMETL
        SSRTMETLV
        SRTMETLVP


3. **gene_expression_detail.gz**: File containing gene expression information for each gene and sample. The output file is a table containing the coding genes in the rows and the samples in the columns. The gene expresison is displayed for each combination. An output example for this tutorial is:

    * Output example:*

    +-------+-------------------------+-----------------------+-------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+--------------------------+-----------------------+------------------------+------------------------+-----------------------+------------------------+-------------------------+------------------------+------------------------+-------------------------+-----------------------+
    | gene  |    simulatedIpp1sample1 |   simulatedIpp1sample2|    simulatedIpp1sample3 |   simulatedIpp1sample4 |   simulatedIpp1sample5 |   simulatedIpp1sample6 |   simulatedIpp1sample7 |   simulatedIpp1sample8 |   simulatedIpp1sample9 |    simulatedIpp1sample10 |  simulatedIpp1sample11|   simulatedIpp1sample12|   simulatedIpp1sample13|  simulatedIpp1sample14|   simulatedIpp1sample15|   simulatedIpp1sample16 |  simulatedIpp1sample17 |  simulatedIpp1sample18 |   simulatedIpp1sample19 |  simulatedIpp1sample20|
    +=======+=========================+=======================+=========================+========================+========================+========================+========================+========================+========================+==========================+=======================+========================+========================+=======================+========================+=========================+========================+========================+=========================+=======================+
    | gene1 |  20688.0                | 17791.0               |         33285.0         | 23488.0                | 46584.0                | 31986.0                | 20888.0                | 32585.0                | 5499.0                 |      13193.0             | 13595.0               |    11495.0             |   3199.0               |   25493.0             |  1899.0                | 10395.0                 |8496.0                  | 20993.0                | 22495.0                 |5199.0                 |
    +-------+-------------------------+-----------------------+-------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+--------------------------+-----------------------+------------------------+------------------------+-----------------------+------------------------+-------------------------+------------------------+------------------------+-------------------------+-----------------------+
    | gene3 |  29585.0                | 58264.0               |         56063.0         | 36978.0                | 85934.0                | 52469.0                | 28083.0                | 51664.0                | 12189.0                |      21579.0             | 33179.0               |    32075.0             |   3496.0               |   31883.0             |  7594.0                | 31576.0                 |16291.0                 | 54973.0                | 23487.0                 |9587.0                 |
    +-------+-------------------------+-----------------------+-------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+--------------------------+-----------------------+------------------------+------------------------+-----------------------+------------------------+-------------------------+------------------------+------------------------+-------------------------+-----------------------+
    | gene4 |  13643.0                | 28684.0               |         57713.0         | 46866.0                | 34632.0                | 43720.0                | 22390.0                | 44068.0                | 11195.0                |      16097.0             | 27289.0               |    17843.0             |   5592.0               |   19243.0             |  9098.0                | 51766.0                 |33586.0                 | 12245.0                |  5948.0                 |19243.0                |
    +-------+-------------------------+-----------------------+-------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+--------------------------+-----------------------+------------------------+------------------------+-----------------------+------------------------+-------------------------+------------------------+------------------------+-------------------------+-----------------------+

4. **ref_sample_peptides.fa.gz**: Fasta file containing the foreground peptides obtained by traversing the splice graph and identifying "short-range" novelty. This file is obtained as output because the command *--output-fasta* is passed to the program. The file contains a header that is the transcript id and the sequence of the corresponding peptide. The name also shows the mutation mode, which in this case is reference. *Note:* As this genome is simulated, there is a higher frequency of stop codons than in nature, that explains the existence of some short peptides.

    *Output example*

    .. code-block::

        fasta
        >gene3:1_3:0:1302:2-exons
        VSDGWACRGSATARPPNPRRAVLCKSIEPTYGRPSV
        >gene7:3_1:0:4980:2-exons
        FGRVPC
        >gene10:1_6:0:6112:2-exons
        YPYVIASCQREAGCRGIICS
        >gene15:0_2:0:9061:2-exons
        LS
        >gene10:1_5:0:6111:2-exons
        ISLCDRKLSEGGGLSRYNLLINCRKGFLGVINRTHVHSLPFRVLIILEPATSLDFRQPGTIDARHCFTMLTGIGNRG
        >gene10:0_7:0:6061:2-exons
        LSLVHPGTRRITKRRRQYPYVIASCQREAGCRGIICS

    This example contains several important things to note.

        - First of all, it is important understand the transcript id. More information about it can be obtained in the :ref:`metadata output file <output-10-build>` section. In this example, the variant id will always be 0 because there is no mutation information.
        - Secondly, it is important to note some short peptides such as *gene7:3_1:0:4980:2-exons* or *gene15:0_2:0:9061:2-exons*. The corresponding peptides are shorter than 9 amino acids, so they will not be shown in the kmer file. This is happening because the translation encountered a stop codon.
        - Finally, in this example, we get three peptides coming from gene 10. However, they are made from different vertex combinations, which results in different peptide sequences.

5. **ref_graph_kmer_SegmExpr:** This is a folder with different files. Each file contains information for the kmers derived from a specific gene, as well as the expression levels of the kmers in each sample. Kmers shown in this file are generated from a single exon and are not located across an exon junction. For a detailed description of the different fields of this file one can refer to the :ref:`file 5 of the output section <output-5-build>`.

    *Output example:*

       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | kmer      | coord                           | isCrossJunction | junctionAnnotated | readFrameAnnotated | simulatedIpp1sample1 | simulatedIpp1sample2 | simulatedIpp1sample3 | simulatedIpp1sample4 | simulatedIpp1sample5 | simulatedIpp1sample6 | simulatedIpp1sample7 | simulatedIpp1sample8 | simulatedIpp1sample9 | simulatedIpp1sample10 | simulatedIpp1sample11 | simulatedIpp1sample12 | simulatedIpp1sample13 | simulatedIpp1sample14 | simulatedIpp1sample15 | simulatedIpp1sample16 | simulatedIpp1sample17 | simulatedIpp1sample18 | simulatedIpp1sample19 | simulatedIpp1sample20 |
       +===========+=================================+=================+===================+====================+======================+======================+======================+======================+======================+======================+======================+======================+======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+
       | RTMETLVP  | 5067:5094:nan:nan:None:None     | False           | False             | True               | 45.64                | 115.36               | 105.56               | 47.54                | 91.91                | 101.54               | 61.33                | 54.53                | 18.4                 | 44.64                 | 72.92                 | 96.13                 | 5.88                  | 90.36                 | 6.78                  | 43.46                 | 20.87                 | 95.22                 |  54.16                | 22.57                 |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | TSSRTMETL | 5061:5088:nan:nan:None:None     | False           | False             | True               | 45.64                | 115.36               | 105.56               | 47.54                | 91.91                | 101.54               | 61.33                | 54.53                | 18.4                 | 44.64                 | 72.92                 | 96.13                 | 5.88                  | 90.36                 | 6.78                  | 43.46                 | 20.87                 | 95.22                 |  54.16                | 22.57                 |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | SSRTMETLV | 5064:5091:nan:nan:None:None     | False           | False             | True               | 45.64                | 115.36               | 105.56               | 47.54                | 91.91                | 101.54               | 61.33                | 54.53                | 18.4                 | 44.64                 | 72.92                 | 96.13                 | 5.88                  | 90.36                 | 6.78                  | 43.46                 | 20.87                 | 95.22                 |  54.16                | 22.57                 |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | LFSDAIRTS | 4063:4090:nan:nan:None:None     | False           | False             | True               | 56.17                | 101.9                | 112.0                | 127.29               | 142.99               | 112.17               | 53.95                | 103.04               | 35.6                 | 54.26                 | 59.72                 | 83.53                 | 10.57                 | 86.79                 | 11.46                 | 51.51                 | 37.02                 |  133.1                |   73.39               |   28.55               |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | HLFSDAIRT | 4060:4087:nan:nan:None:None     | False           | False             | True               | 56.17                | 101.9                | 112.0                | 127.29               | 142.99               | 112.17               | 53.95                | 103.04               | 35.6                 | 54.26                 | 59.72                 | 83.53                 | 10.57                 | 86.79                 | 11.46                 | 51.51                 | 37.02                 | 133.1                 |   73.39               |   28.55               |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+

   This example contains several important things to note.

        - The folder contains 9 files. This is because in the file *ref_sample_peptides.fa.gz* there are peptides for the 9 different genes.
        - In this example the results for two different genes are shown. One can see that by looking at the expression levels. The expression levels for the first three kmers are the same, and the same happens for the two last kmers. This is because the gene expression is obtained at a per-gene basis. Therefore, all kmers derived from the same gene will have the same expression level.
        - As we are dealing with segment kmers, the fields isCrossJunction and junctionAnnotated are always False.
        - The field readFrameAnnotated shows whether the kmers were obtained from a read frame present in the annotation file or if they were obtained by reading frame propagation.

6. **ref_graph_kmer_JuncExpr:** This is a folder containing different files. In this case, it there are three files. Each file shows the expression levels for different kmers, across the 20 samples. For a detailed description of the different fields of this file one can refer to the :ref:`file 6 of the output section <output-6-build>`.

    *Output example:*

       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | kmer      | coord                           | isCrossJunction | junctionAnnotated | readFrameAnnotated | simulatedIpp1sample1 | simulatedIpp1sample2 | simulatedIpp1sample3 | simulatedIpp1sample4 | simulatedIpp1sample5 | simulatedIpp1sample6 | simulatedIpp1sample7 | simulatedIpp1sample8 | simulatedIpp1sample9 | simulatedIpp1sample10 | simulatedIpp1sample11 | simulatedIpp1sample12 | simulatedIpp1sample13 | simulatedIpp1sample14 | simulatedIpp1sample15 | simulatedIpp1sample16 | simulatedIpp1sample17 | simulatedIpp1sample18 | simulatedIpp1sample19 | simulatedIpp1sample20 |
       +===========+=================================+=================+===================+====================+======================+======================+======================+======================+======================+======================+======================+======================+======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+
       | SSSLVSDGW | 1140:1150:1300:1317:None:None   |   True          |  False            |      True          |    92.0              |  183.0               |  175.0               |  119.0               | 285.0                | 153.0                | 85.0                 | 170.0                |  40.0                |  54.0                 |  95.0                 |  99.0                 | 10.0                  |  93.0                 |  24.0                 |  87.0                 |  41.0                 |  172.0                |  69.0                 | 32.0                  |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | EPPTYGRPSV| 1383:1400:1500:1510:None:None   |   True          |  False            |      False         |    54.0              |  117.0               |  75.0                |  77.0                | 138.0                | 70.0                 | 33.0                 | 93.0                 |  23.0                |  31.0                 |  46.0                 |  75.0                 | 6.0                   |  42.0                 |  4.0                  |  61.0                 |  26.0                 |  83.0                 |  30.0                 | 6.0                   |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       |RESSSLVSD  | 1134:1150:1300:1311:None:None   |   True          |  False            |      True          |    92.0              |  183.0               |  175.0               |  119.0               | 285.0                | 153.0                | 85.0                 | 170.0                |  40.0                |  54.0                 |  95.0                 |  99.0                 | 10.0                  |  93.0                 |  24.0                 |  87.0                 |  41.0                 |  172.0                |  69.0                 | 32.0                  |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+


In this file, all the kmers appearing belong to a junction between exons. Therefore, the field *isCrossJunction* will always have the value True.

7. **expression_counts.libsize.tsv:** File containing 75% of expression and total expression for each sample. For a given sample, the “75% of expression” is defined as the 75th quantile of the gene expression distribution across coding genes. For a given sample the “total expression” is defined as the total gene expression across coding genes. Generated only if *–-disable-libsize* is set to False and if *-–count-path* file is provided. It is computed from the file *gene_expression_detail.gz*.

    *Output example:*

    +-----------------------+-------------------------+------------------------+
    |         sample        |    libsize_75percent    |   libsize_total_count  |
    +-----------------------+-------------------------+------------------------+
    |  simulatedIpp1sample1 |         29374.0         |        212925.0        |
    +-----------------------+-------------------------+------------------------+
    |  simulatedIpp1sample2 |         43077.0         |        361533.0        |
    +-----------------------+-------------------------+------------------------+
    |  simulatedIpp1sample3 |         57713.0         |        485196.0        |
    +-----------------------+-------------------------+------------------------+
    |  simulatedIpp1sample4 |         51055.0         |        383447.0        |
    +-----------------------+-------------------------+------------------------+
    |  simulatedIpp1sample5 |         73240.0         |        541376.0        |
    +-----------------------+-------------------------+------------------------+
    |  simulatedIpp1sample6 |         52469.0         |        417486.0        |
    +-----------------------+-------------------------+------------------------+
    |  simulatedIpp1sample7 |         26440.0         |        239553.0        |
    +-----------------------+-------------------------+------------------------+
    |  simulatedIpp1sample8 |         44205.0         |        390951.0        |
    +-----------------------+-------------------------+------------------------+
    |  simulatedIpp1sample9 |         13190.0         |        103124.0        |
    +-----------------------+-------------------------+------------------------+
    | simulatedIpp1sample10 |         20537.0         |        180067.0        |
    +-----------------------+-------------------------+------------------------+
    | simulatedIpp1sample11 |         45020.0         |        305882.0        |
    +-----------------------+-------------------------+------------------------+
    | simulatedIpp1sample12 |         37876.0         |        300323.0        |
    +-----------------------+-------------------------+------------------------+
    | simulatedIpp1sample13 |         6248.0          |         48668.0        |
    +-----------------------+-------------------------+------------------------+
    | simulatedIpp1sample14 |         40640.0         |        339159.0        |
    +-----------------------+-------------------------+------------------------+
    | simulatedIpp1sample15 |         6696.0          |         49765.0        |
    +-----------------------+-------------------------+------------------------+
    | simulatedIpp1sample16 |         24288.0         |        196125.0        |
    +-----------------------+-------------------------+------------------------+
    | simulatedIpp1sample17 |         12985.0         |        110283.0        |
    +-----------------------+-------------------------+------------------------+
    | simulatedIpp1sample18 |         54973.0         |        442504.0        |
    +-----------------------+-------------------------+------------------------+
    | simulatedIpp1sample19 |         33586.0         |        263805.0        |
    +-----------------------+-------------------------+------------------------+
    | simulatedIpp1sample20 |         13190.0         |        109355.0        |
    +-----------------------+-------------------------+------------------------+

8. **Annot_IS_SUCCESS:** This is an empty file. It is obtained because the generation of the background (or annotation) files was successful. If the generation of the background files was not successful, this file would not be generated.

9. **output_sample_IS_SUCCESS:** This is an empty file. It is obtained because the generation of the foreground (or sample) files was successful. If the generation of the foreground files was not successful, this file would not be generated.

10. **somatic_and_germline_sample_peptides_meta.gz:** File containing details for each peptide. A detailed explanation of the output can be seen in :ref:`metadata output file <output-10-build>`.

    *Output example:*

    +--------------------------------------------------------------------------------------------+-------------------------------+---------------------+---------------------+----------+------------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+
    |                                peptide                                                     |              id               | readFrame           | readFrameAnnotated  | geneName | geneChr    | geneStrand |     mutationMode      | hasStopCodon | isInJunctionList  | isIsolated | variantComb  | variantSegExpr  | modifiedExonsCoord  | originalExonsCoord  | vertexIdx |   kmerType   |
    +--------------------------------------------------------------------------------------------+-------------------------------+---------------------+---------------------+----------+------------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+
    | VSDGWACRGSATARPPNPRRAVLCKSIEPTYGRPSV                                                       | gene3:1_3:0:1302:2-exons      | 2                   | False               | gene3    | 1          | +          | ref                   | 1            | nan               | 0          | nan          | nan             | 1302;1400;1500;1600 | 1300;1400;1500;1600 | 1;3       | 2-exons      |
    +--------------------------------------------------------------------------------------------+-------------------------------+---------------------+---------------------+----------+------------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+
    | FGRVPC                                                                                     | gene7:3_1:0:4980:2-exons      | 2                   | True                | gene7    | 1          | -          | ref                   | 1            | nan               | 1          | nan          | nan             | 4900;4980;4700;4800 | 4900;5000;4700;4800 | 3;1       | 2-exons      |
    +--------------------------------------------------------------------------------------------+-------------------------------+---------------------+---------------------+----------+------------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+
    | YPYVIASCQREAGCRGIICS                                                                       | gene10:1_6:0:6112:2-exons     | 0                   | True                | gene10   | 1          | +          | ref                   | 1            | nan               | 1          | nan          | nan             | 6112;6250;6400;6748 | 6100;6250;6400;6750 | 1;6       | 2-exons      |
    +--------------------------------------------------------------------------------------------+-------------------------------+---------------------+---------------------+----------+------------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+
    | YPYVIASCQREAGCRGIICS                                                                       | gene10:1_5:0:6112:2-exons     | 0                   | True                | gene10   | 1          | +          | ref                   | 1            | nan               | 1          | nan          | nan             | 6112;6250;6400;6498 | 6100;6250;6400;6500 | 1;5       | 2-exons      |
    +--------------------------------------------------------------------------------------------+-------------------------------+---------------------+---------------------+----------+------------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+
    | ISLCDRKLSEGGGLSRYNLLINCRKGFLGVINRTHVHSLPFRVLIILEPATSLDFRQPGTIDARHCFTMLTGIGNRG              | gene10:1_5:0:6111:2-exons     | 1                   | True                | gene10   | 1          | +          | ref                   | 1            | nan               | 0          | nan          | nan             | 6111;6250;6400;6498 | 6100;6250;6400;6500 | 1;5       | 2-exons      |
    +--------------------------------------------------------------------------------------------+-------------------------------+---------------------+---------------------+----------+------------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+

    Things to note from the table above:

        - All the ids have a variant number equal to 0. The same happens with VariantComb and VariantSegExpr, which are nan. This is because mutations are not provided.
        - In the third and fourth column, one can see two peptides that have the same sequence but come from different vertices. Moreover, they both have the same reading frame. This is because the sequence has a stop codon in exon 1 if reading frame 0 is used, so that in both cases only the part of exon 1 up to the mutation is translated.
        - On the other hand, on the fifth line, we can see how a shift of 1 nucleotide in the reading frame leads to the disappearance of the stop codon, so that the whole sequence is translated.


Build mode *with* mutations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the second part of the example, we introduce somatic and germline mutations in the analysis. The command used in this tutorial to run the *build* mode is:

.. code-block:: console

    immunopepper  build --output-dir immunopepper_usecase/ --ann-path  immunopepper/tests/data_simulated/data/simulated_Ipp.gtf --splice-path  immunopepper/tests/data_simulated/data/genes_graph_conf3.merge_graphs.pickle --ref-path  immunopepper/tests/data_simulated/data/genome.fa --kmer 9 --count-path immunopepper/tests/data_simulated/data/genes_graph_conf3.merge_graphs.count.hdf5 --parallel 1 --batch-size 1  --start-id 0 --process-num 0 --output-fasta --somatic immunopepper/tests/data_simulated/data/variants_somatic.vcf --germline immunopepper/tests/data_simulated/data/variants_germline.vcf --mutation-sample simulated_Ipp_1_sample3 --verbose 2

#TODO: save the results in the github so that users can look at them without running the example?

In this command, the build mode of immunopepper is run on the :ref:`input_data` described in the section above. Moreover, the output directory is set to a folder called *immunopepper_usecase/cohort_mutsimulated_Ipp_1_sample3*, located on the directory where the command is executed. The kmer length is set to 9, as it is a common kmer length selected in clinical applications. Finally, there are also two mutation files provided, a somatic and a germline file. These files will apply the existing mutations and take them into account when computing the output.
One important thing to note is that, if mutations are provided, an extra filter layer is included. This layer will ensure that only peptides different to the reference (base genome + germline) are included in the output.

**Terminal output:**

The output displayed in the terminal is the following:

.. code-block:: console

    2023-06-20 19:21:51,580 INFO     Command lineNamespace(output_dir='immunopepper_usecase/', ann_path='immunopepper/tests/data_simulated/data/simulated_Ipp.gtf', splice_path='immunopepper/tests/data_simulated/data/genes_graph_conf3.merge_graphs.pickle', ref_path='immunopepper/tests/data_simulated/data/genome.fa', kmer=9, libsize_extract=False, all_read_frames=False, count_path='immunopepper/tests/data_simulated/data/genes_graph_conf3.merge_graphs.count.hdf5', output_samples=[], heter_code=0, compressed=True, parallel=1, batch_size=1, pickle_samples=[], process_chr=None, complexity_cap=None, genes_interest=None, start_id=0, process_num=0, skip_annotation=False, keep_tmpfiles=False, libsize_path=None, output_fasta=True, force_ref_peptides=False, filter_redundant=False, kmer_database=None, gtex_junction_path=None, disable_concat=False, disable_process_libsize=False, mutation_sample='simulated_Ipp_1_sample3', germline='immunopepper/tests/data_simulated/data/variants_germline.vcf', somatic='immunopepper/tests/data_simulated/data/variants_somatic.vcf', sample_name_map=None, use_mut_pickle=False, verbose=2)
    2023-06-20 19:21:51,580 INFO     >>>>>>>>> Build: Start Preprocessing
    2023-06-20 19:21:51,580 INFO     Building lookup structure ...
    2023-06-20 19:21:51,581 INFO            Time spent: 0.000 seconds
    2023-06-20 19:21:51,581 INFO            Memory usage: 0.146 GB
    2023-06-20 19:21:51,581 INFO     Loading count data ...
    2023-06-20 19:21:51,584 INFO            Time spent: 0.003 seconds
    2023-06-20 19:21:51,584 INFO            Memory usage: 0.147 GB
    2023-06-20 19:21:51,585 INFO     Loading splice graph ...
    2023-06-20 19:21:51,586 INFO            Time spent: 0.000 seconds
    2023-06-20 19:21:51,586 INFO            Memory usage: 0.148 GB
    2023-06-20 19:21:51,586 INFO     Add reading frame to splicegraph ...
    2023-06-20 19:21:51,588 INFO            Time spent: 0.002 seconds
    2023-06-20 19:21:51,588 INFO            Memory usage: 0.148 GB
    2023-06-20 19:21:51,588 INFO     >>>>>>>>> Finish Preprocessing
    2023-06-20 19:21:51,588 INFO     >>>>>>>>> Start traversing splicegraph
    2023-06-20 19:21:51,588 INFO     >>>> Processing output_sample cohort, there are 9 graphs in total
    2023-06-20 19:21:51,588 INFO     Saving results to immunopepper_usecase/cohort_mutsimulated_Ipp_1_sample3
    2023-06-20 19:21:51,588 INFO     Not Parallel
    2023-06-20 19:21:51,588 INFO     >>>>>>>>> Start Background processing
    2023-06-20 19:21:51,591 INFO     Saved somatic_and_germline_annot_peptides.fa.gz with 40 lines in 0.0004s
    2023-06-20 19:21:51,591 INFO     Saved somatic_and_germline_annot_kmer.gz with 298 lines in 0.0003s
    2023-06-20 19:21:51,592 DEBUG    ....cohort: annotation graph from batch all/9 processed, max time cost: 0.0, memory cost: 0.15 GB
    2023-06-20 19:21:51,592 INFO     >>>>>>>>> Start Foreground processing
    2023-06-20 19:21:51,632 INFO     Saved gene_expression_detail.gz with 9 lines in 0.0004s
    2023-06-20 19:21:51,632 INFO     Saved somatic_and_germline_sample_peptides.fa.gz with 46 lines in 0.0003s
    2023-06-20 19:21:51,633 INFO     Saved somatic_and_germline_sample_peptides_meta.gz with 23 lines in 0.0004s
    2023-06-20 19:21:51,633 DEBUG    ....cohort: output_sample graph from batch all/9 processed, max time cost: 0.01, memory cost: 0.15 GB
    2023-06-20 19:21:51,645 INFO     Saved library size results to immunopepper_usecase/expression_counts.libsize.tsv

**Output files:**


1. **somatic_and_germline_annot_peptides.fa.gz**: This is a fasta file containing the background peptides for each gene transcript. The file contains a header that is the transcript id and the sequence of the corresponding peptide. The name also shows the mutation mode, which in this case is somatic and germline. *Note*: As this genome is simulated, there is a higher frequency of stop codons than in nature, that explains the existence of some short peptides.

    *Output example:*

    .. code-block::

        >gene8.t1
        TSSRTMETLVP
        >gene4.t2
        SVRRTPRFRRTAEAPVSRSLIITHLGDGGWEP
        >gene15.t1
        TFLKKSLRLNSI

    The peptides of the annotation will aready have the germline variants included.

2. **somatic_and_germline_annot_kmer.gz**: This file contains the kmers of length 9 obtained from the background peptides. The kmers are obtained by passing a sliding window through the peptides contained in somatic_and_germline_annot_peptides.fa.gz. *Note*: For that peptides that are shorter than 9, the kmers are not obtained.

    The kmers in this example are obtained from *gene8.t1* and *gene4.t2*.

    *Output example:*

    .. code-block::

        kmer
        TSSRTMETL
        SSRTMETLV
        SRTMETLVP
        SVRRTPRFR
        VRRTPRFRRT


3. **gene_expression_detail.gz**: File containing gene expression information for each gene and sample. The output file is a table containing the coding genes in the rows and the samples in the columns. The gene expresison is displayed for each combination. An output example for this tutorial is:

    * Output example:*

    +-------+-------------------------+-----------------------+-------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+--------------------------+-----------------------+------------------------+------------------------+-----------------------+------------------------+-------------------------+------------------------+------------------------+-------------------------+-----------------------+
    | gene  |    simulatedIpp1sample1 |   simulatedIpp1sample2|    simulatedIpp1sample3 |   simulatedIpp1sample4 |   simulatedIpp1sample5 |   simulatedIpp1sample6 |   simulatedIpp1sample7 |   simulatedIpp1sample8 |   simulatedIpp1sample9 |    simulatedIpp1sample10 |  simulatedIpp1sample11|   simulatedIpp1sample12|   simulatedIpp1sample13|  simulatedIpp1sample14|   simulatedIpp1sample15|   simulatedIpp1sample16 |  simulatedIpp1sample17 |  simulatedIpp1sample18 |   simulatedIpp1sample19 |  simulatedIpp1sample20|
    +=======+=========================+=======================+=========================+========================+========================+========================+========================+========================+========================+==========================+=======================+========================+========================+=======================+========================+=========================+========================+========================+=========================+=======================+
    | gene1 |  20688.0                | 17791.0               |         33285.0         | 23488.0                | 46584.0                | 31986.0                | 20888.0                | 32585.0                | 5499.0                 |      13193.0             | 13595.0               |    11495.0             |   3199.0               |   25493.0             |  1899.0                | 10395.0                 |8496.0                  | 20993.0                | 22495.0                 |5199.0                 |
    +-------+-------------------------+-----------------------+-------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+--------------------------+-----------------------+------------------------+------------------------+-----------------------+------------------------+-------------------------+------------------------+------------------------+-------------------------+-----------------------+
    | gene3 |  29585.0                | 58264.0               |         56063.0         | 36978.0                | 85934.0                | 52469.0                | 28083.0                | 51664.0                | 12189.0                |      21579.0             | 33179.0               |    32075.0             |   3496.0               |   31883.0             |  7594.0                | 31576.0                 |16291.0                 | 54973.0                | 23487.0                 |9587.0                 |
    +-------+-------------------------+-----------------------+-------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+--------------------------+-----------------------+------------------------+------------------------+-----------------------+------------------------+-------------------------+------------------------+------------------------+-------------------------+-----------------------+
    | gene4 |  13643.0                | 28684.0               |         57713.0         | 46866.0                | 34632.0                | 43720.0                | 22390.0                | 44068.0                | 11195.0                |      16097.0             | 27289.0               |    17843.0             |   5592.0               |   19243.0             |  9098.0                | 51766.0                 |33586.0                 | 12245.0                |  5948.0                 |19243.0                |
    +-------+-------------------------+-----------------------+-------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+--------------------------+-----------------------+------------------------+------------------------+-----------------------+------------------------+-------------------------+------------------------+------------------------+-------------------------+-----------------------+

4. **somatic_and_germline_sample_peptides.fa.gz**: Fasta file containing the foreground peptides obtained by traversing the splice graph and taking into account the mutations. This file is obtained as output because the command *--output-fasta* is passed to the program. The file contains a header that is the transcript id and the sequence of the corresponding peptide. The name also shows the mutation mode, which in this case is somatic and germline. *Note:* As this genome is simulated, there is a higher frequency of stop codons than in nature, that explains the existence of some short peptides.

    *Output example*

    .. code-block::

        >gene4:0_1:3:2112:2-exons
        SVRRTPRFRRTAEAPVSRSLIITHLGDEGWEP
        >gene4:0_2:6:2112:2-exons
        SVRRTPRFRRTAEAPVSRSLIITHLGDEGWEP
        >gene4:0_1:5:2112:2-exons
        SVRRTPRFRRTAEAPVSRSLIITHLGDEGWEP
        >gene3:0_1:0:1062:2-exons
        NEVIGECIACSASFDATTIGRSRHRESSSLVSDGWACRGSATARPPNPRRAVLCKSIEPTYA
        >gene3:0_2:3:1062:2-exons
        NEVIGECIACSASFDATTIGRSRHRESSSLVSDGWACRGSATARPPNPRRAVLCKSIEPTYAR

    Important things to note:
        - When the somatic and germline files are included, only peptides belonging to gene4 and gene3 are given as output. This is because only the peptides that are different to the reference are given as an output. For many of the genes and exons, there are not mutations present, which means that the peptides will be the same as in the reference.
        - The peptides coming from gene4 are equal to each other. This is because there is a stop codon in exon 0, and the peptide is truncated there. It can be checked by looking at the isIsolated field in the metadata.

5. **somatic_and_germline_graph_kmer_SegmExpr:** This is a folder with different files. Each file contains information for different kmers derived from *somatic_and_germline_sample_peptides.fa.gz*. The files contain information about the expression levels of kmers found in an exon. Kmers shown in this file are generated from a single exon and are not located across an exon junction. For a detailed description of the different fields of this file one can refer to the :ref:`file 5 of the output section <output-5-build>`.

    In this example one can see three peptides derived from the same gene. This explains why the expression level is the same in each sample for the three kmers, as gene expression is computed per gene.

    *Output example:*

       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | kmer      | coord                           | isCrossJunction | junctionAnnotated | readFrameAnnotated | simulatedIpp1sample1 | simulatedIpp1sample2 | simulatedIpp1sample3 | simulatedIpp1sample4 | simulatedIpp1sample5 | simulatedIpp1sample6 | simulatedIpp1sample7 | simulatedIpp1sample8 | simulatedIpp1sample9 | simulatedIpp1sample10 | simulatedIpp1sample11 | simulatedIpp1sample12 | simulatedIpp1sample13 | simulatedIpp1sample14 | simulatedIpp1sample15 | simulatedIpp1sample16 | simulatedIpp1sample17 | simulatedIpp1sample18 | simulatedIpp1sample19 | simulatedIpp1sample20 |
       +===========+=================================+=================+===================+====================+======================+======================+======================+======================+======================+======================+======================+======================+======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+
       | RFRRTAEAP |    2130:2157:nan:nan:None:None  |   False         | False             |  True              |  36.98               |  79.75               | 161.4                | 128.52               | 101.0                |  119.43              | 62.57                | 119.59               | 29.76                | 48.62                 |  77.01                | 51.13                 |  15.85                |  125.49               | 15.11                 | 54.82                 | 23.91                 | 141.25                | 95.49                 | 33.13                 |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | APVSRSLII |    2151:2178:nan:nan:None:None  |   False         | False             |  True              |  36.98               |  79.75               | 161.4                | 128.52               | 101.0                |  119.43              | 62.57                | 119.59               | 29.76                | 48.62                 |  77.01                | 51.13                 |  15.85                |  125.49               | 15.11                 | 54.82                 | 23.91                 | 141.25                | 95.49                 | 33.13                 |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | RRTPRFRRT |    2118:2145:nan:nan:None:None  |   False         | False             |  True              |  36.98               |  79.75               | 161.4                | 128.52               | 101.0                |  119.43              | 62.57                | 119.59               | 29.76                | 48.62                 |  77.01                | 51.13                 |  15.85                |  125.49               | 15.11                 | 54.82                 | 23.91                 | 141.25                | 95.49                 | 33.13                 |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+



6. **somatic_and_germline_graph_kmer_JuncExpr:** This is a folder containing different files. Each file has information for different kmers derived from *somatic_and_germline_sample_peptides.fa.gz*. The files contain information about the expression levels of kmers found across an exon junction. For a detailed description of the different fields of this file one can refer to the :ref:`file 6 of the output section <output-6-build>`.

    *Output example:*

       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | kmer      | coord                           | isCrossJunction | junctionAnnotated | readFrameAnnotated | simulatedIpp1sample1 | simulatedIpp1sample2 | simulatedIpp1sample3 | simulatedIpp1sample4 | simulatedIpp1sample5 | simulatedIpp1sample6 | simulatedIpp1sample7 | simulatedIpp1sample8 | simulatedIpp1sample9 | simulatedIpp1sample10 | simulatedIpp1sample11 | simulatedIpp1sample12 | simulatedIpp1sample13 | simulatedIpp1sample14 | simulatedIpp1sample15 | simulatedIpp1sample16 | simulatedIpp1sample17 | simulatedIpp1sample18 | simulatedIpp1sample19 | simulatedIpp1sample20 |
       +===========+=================================+=================+===================+====================+======================+======================+======================+======================+======================+======================+======================+======================+======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+
       | RHRESSSLV |   1128:1150:1300:1305:None:None |   True          | False             |  True              |  92.0                |  183.0               | 175.0                | 119.0                | 285.0                | 153.0                | 85.0                 | 170.0                | 40.0                 | 54.0                  | 95.0                  | 99.0                  | 10.0                  | 93.0                  | 24.0                  | 87.0                  | 41.0                  | 172.0                 | 69.0                  | 32.0                  |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | EPTYARPSV |   1383:1400:1500:1510:None:None |   True          | False             |  True              |  54.0                |  117.0               | 75.0                 | 77.0                 | 138.0                | 70.0                 | 33.0                 | 93.0                 | 23.0                 | 31.0                  | 46.0                  | 75.0                  | 6.0                   | 42.0                  | 4.0                   | 61.0                  | 26.0                  | 83.0                  | 30.0                  | 6.0                   |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+
       | KSIEPTYAR |   1374:1400:1500:1501:None:None |  True           | False             |  False             |  54.0                |  117.0               | 75.0                 | 77.0                 | 138.0                | 70.0                 | 33.0                 | 93.0                 | 23.0                 | 31.0                  | 46.0                  | 75.0                  | 6.0                   | 42.0                  | 4.0                   | 61.0                  | 26.0                  | 83.0                  | 30.0                  | 6.0                   |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+


7. **expression_counts.libsize.tsv:** File containing 75% of expression and total expression for each sample. For a given sample, the “75% of expression” is defined as the 75th quantile of the gene expression distribution across coding genes. For a given sample the “total expression” is defined as the total gene expression across coding genes. Generated only if –disable-libsize is set to False and if –count-path file is provided. It is computed from the file gene_expression_detail.gz.

    The output is the same as in the tutorial without mutations. This is because this file,as well as the file gene_expression_detail.gz, are generated from the same count file. Therefore, including mutations does not have an influence on the output of this file.

    *Output example:*

        +-----------------------+-------------------------+------------------------+
        |         sample        |    libsize_75percent    |   libsize_total_count  |
        +-----------------------+-------------------------+------------------------+
        |  simulatedIpp1sample1 |         29374.0         |        212925.0        |
        +-----------------------+-------------------------+------------------------+
        |  simulatedIpp1sample2 |         43077.0         |        361533.0        |
        +-----------------------+-------------------------+------------------------+
        |  simulatedIpp1sample3 |         57713.0         |        485196.0        |
        +-----------------------+-------------------------+------------------------+
        |  simulatedIpp1sample4 |         51055.0         |        383447.0        |
        +-----------------------+-------------------------+------------------------+
        |  simulatedIpp1sample5 |         73240.0         |        541376.0        |
        +-----------------------+-------------------------+------------------------+
        |  simulatedIpp1sample6 |         52469.0         |        417486.0        |
        +-----------------------+-------------------------+------------------------+
        |  simulatedIpp1sample7 |         26440.0         |        239553.0        |
        +-----------------------+-------------------------+------------------------+
        |  simulatedIpp1sample8 |         44205.0         |        390951.0        |
        +-----------------------+-------------------------+------------------------+
        |  simulatedIpp1sample9 |         13190.0         |        103124.0        |
        +-----------------------+-------------------------+------------------------+
        | simulatedIpp1sample10 |         20537.0         |        180067.0        |
        +-----------------------+-------------------------+------------------------+
        | simulatedIpp1sample11 |         45020.0         |        305882.0        |
        +-----------------------+-------------------------+------------------------+
        | simulatedIpp1sample12 |         37876.0         |        300323.0        |
        +-----------------------+-------------------------+------------------------+
        | simulatedIpp1sample13 |         6248.0          |         48668.0        |
        +-----------------------+-------------------------+------------------------+
        | simulatedIpp1sample14 |         40640.0         |        339159.0        |
        +-----------------------+-------------------------+------------------------+
        | simulatedIpp1sample15 |         6696.0          |         49765.0        |
        +-----------------------+-------------------------+------------------------+
        | simulatedIpp1sample16 |         24288.0         |        196125.0        |
        +-----------------------+-------------------------+------------------------+
        | simulatedIpp1sample17 |         12985.0         |        110283.0        |
        +-----------------------+-------------------------+------------------------+
        | simulatedIpp1sample18 |         54973.0         |        442504.0        |
        +-----------------------+-------------------------+------------------------+
        | simulatedIpp1sample19 |         33586.0         |        263805.0        |
        +-----------------------+-------------------------+------------------------+
        | simulatedIpp1sample20 |         13190.0         |        109355.0        |
        +-----------------------+-------------------------+------------------------+

8. **Annot_IS_SUCCESS:** This is an empty file. It is obtained because the generation of the background (or annotation) files was successful. If the generation of the background files was not successful, this file would not be generated.

9. **output_sample_IS_SUCCESS:** This is an empty file. It is obtained because the generation of the foreground (or sample) files was successful. If the generation of the foreground files was not successful, this file would not be generated.

10. **somatic_and_germline_sample_peptides_meta.gz:** File containing details for each peptide. A detailed explanation of the output can be seen in :ref:`metadata output file <output-10-build>`.

    *Output example:*

    +--------------------------------------------------------------------------------------------+-------------------------------+-----------+-----------------------+----------+----------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+
    |                                peptide                                                     |              id               | readFrame |   readFrameAnnotated  | geneName | geneChr  | geneStrand |     mutationMode      | hasStopCodon | isInJunctionList  | isIsolated | variantComb  | variantSegExpr  | modifiedExonsCoord  | originalExonsCoord  | vertexIdx |   kmerType   |
    +--------------------------------------------------------------------------------------------+-------------------------------+-----------+-----------------------+----------+----------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+
    | SVRRTPRFRRTAEAPVSRSLIITHLGDEGWEP                                                           | gene4:0_1:2:2112:2-exons      |         0 |      True             |   gene4  |    1     |      +     | somatic_and_germline  |       1      |        nan        |     1      |   2194;2249  |       161       | 2112;2250;2350;2449 | 2100;2250;2350;2450 |    0;1    |   2-exons    |
    +--------------------------------------------------------------------------------------------+-------------------------------+-----------+-----------------------+----------+----------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+
    | SVRRTPRFRRTAEAPVSRSLIITHLGDEGWEP                                                           | gene4:0_2:0:2112:2-exons      |         0 |    True               |   gene4  |    1     |      +     | somatic_and_germline  |       1      |        nan        |     1      |     2194     |       161       | 2112;2250;2550;2649 | 2100;2250;2550;2650 |    0;2    |   2-exons    |
    +--------------------------------------------------------------------------------------------+-------------------------------+-----------+-----------------------+----------+----------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+
    | NEVIGECIACSASFDATTIGRSRHRESSSLVSDGWACRGSATARPPNPRRAVLCKSIEPTYAR                            | gene3:0_2:0:1062:2-exons      |         1 |     True              |   gene3  |    1     |      +     | somatic_and_germline  |       1      |        nan        |     0      |     1396     |       177       | 1062;1150;1300;1599 | 1050;1150;1300;1600 |    0;2    |   2-exons    |
    +--------------------------------------------------------------------------------------------+-------------------------------+-----------+-----------------------+----------+----------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+

    Important things to note:

        - In this case, the mutation index in the id is not always 0, since we have a mutation file and we study several mutations.
        - The field *variantComb* shows the variant combination used to create the peptide. Variants in the software are SNP at DNA level.
        - The field *variantSegExpr* shows the expression of the gene with the variant combination used to create the peptide.
        - We can see once more that the sequence for the gene4 peptides is the same. This is because there is a stop codon in the exon 0 of the gene, so that the peptide is truncated.

.. _tutorial_samplespecif:

Running samplespecif mode
--------------------------
