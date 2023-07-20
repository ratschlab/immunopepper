Tutorials
==========

In this page, a list of tutorials is provided to learn how to use the software. First, the input data used for this tutorial will be described. Then, a tutorial for each mode will be presented.

1. :ref:`input_data`: Description of the input data.
2. :ref:`tutorial_build`: Tutorial for the *build* mode.
3. :ref:`tutorial_samplespecif`: Tutorial for the *samplespecif* mode.
4. :ref:`tutorial_cancerspecif`: Tutorial for the *cancerspecif* mode.
5. :ref:`tutorial_mhcbind`: Tutorial for the *mhcbind* mode.

.. note:: The folder with the tutorial results will be generated in the directory where immunopepper is run, not inside the immunopepper folder.

.. _input_data:

Input data
--------------

The input data used for this tutorial can be found in the folder *'tests/data_simulated/data/build_mode'*. The data is composed of the following files:

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

Command
~~~~~~~~~~
First of all, we will show an example with the basic commands, without any additional functionalities such as mutations. The command used for this part of the tutorial is:

.. code-block::

    immunopepper  build --output-dir immunopepper_usecase/ --ann-path  immunopepper/tests/data_simulated/data/build_mode/simulated_Ipp.gtf --splice-path  immunopepper/tests/data_simulated/data/build_mode/genes_graph_conf3.merge_graphs.pickle --ref-path  immunopepper/tests/data_simulated/data/build_mode/genome.fa --kmer 9 --count-path immunopepper/tests/data_simulated/data/build_mode/genes_graph_conf3.merge_graphs.count.hdf5 --parallel 1 --batch-size 1  --start-id 0 --process-num 0 --output-fasta --verbose 2

By calling this command, the software will generate the possible kmers/peptides for each of the 9 genes in the splice graph. It will take into account the reference genome and the annotation file, and it will generate an output for both the background and foreground peptides. The command is run on the :ref:`input_data` described in the section above. Moreover, the output directory is set to a folder called *immunopepper_usecase*, located on the directory where the command is executed. The kmer length is set to 9, as it is a common kmer length selected in clinical applications.

Terminal output:
~~~~~~~~~~~~~~~~~~

The output displayed in the command line is the following:

.. code-block:: console

    2023-06-22 12:48:54,100 INFO     Command lineNamespace(output_dir='immunopepper_usecase/', ann_path='immunopepper/tests/data_simulated/data/build_mode/simulated_Ipp.gtf', splice_path='immunopepper/tests/data_simulated/data/build_mode/genes_graph_conf3.merge_graphs.pickle', ref_path='immunopepper/tests/data_simulated/data/build_mode/genome.fa', kmer=9, libsize_extract=False, all_read_frames=False, count_path='immunopepper/tests/data_simulated/data/build_mode/genes_graph_conf3.merge_graphs.count.hdf5', output_samples=[], heter_code=0, compressed=True, parallel=1, batch_size=1, pickle_samples=[], process_chr=None, complexity_cap=None, genes_interest=None, start_id=0, process_num=0, skip_annotation=False, libsize_path=None, output_fasta=True, force_ref_peptides=False, filter_redundant=False, kmer_database=None, gtex_junction_path=None, disable_concat=False, disable_process_libsize=False, mutation_sample=None, germline='', somatic='', sample_name_map=None, use_mut_pickle=False, verbose=2)
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

Output files
~~~~~~~~~~~~~

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
        ...

2. **ref_annot_kmer.gz**: This file contains the kmers of length 9 obtained from the background peptides. The kmers are obtained by passing a sliding window through the peptides contained in ref_annot_peptides.fa.gz. *Note*: For peptides that are shorter than 9 aminoacids, the kmers are not obtained.

    In the example below, one can see the kmers obtained from *gene8.t1*.

    *Output example:*

    .. code-block::

        kmer
        TSSRTMETL
        SSRTMETLV
        SRTMETLVP
        ...


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
    | ...   |  ...                    | ...                   |         ...             | ...                    | ...                    | ...                    | ...                    | ...                    | ...                    |      ...                 | ...                   |    ...                 |   ...                  |   ...                 |  ...                   | ...                     |...                     | ...                    | ...                     |...                    |
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
        ...

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
       | ...       |  ...                            | ...             |         ...       | ...                | ...                  | ...                  | ...                  | ...                  | ...                  |      ...             | ...                  |    ...               |   ...                |   ...                 |  ...                  | ...                   |...                    | ...                   | ...                   |...                    | ...                   | ...                   | ...                   | ...                   |
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
       | ...       |  ...                            |   ...           |  ...              |       ...          |    ...               |  ...                 |   ...                |  ...                 | ...                  | ...                  | ...                  | ...                  |   ...                |   ...                 |  ...                  | ...                   |...                    | ...                   | ...                   |...                    | ...                   | ...                   | ...                   | ...                   |
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
    | ...                                                                                        | ...                           | ...                 | ...                 | ...      | ...        | ...        | ...                   | ...          | ...               | ...        | ...          | ...             | ...                 | ...                 | ...       | ...          |
    +--------------------------------------------------------------------------------------------+-------------------------------+---------------------+---------------------+----------+------------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+

    Things to note from the table above:

        - All the ids have a variant number equal to 0. The same happens with VariantComb and VariantSegExpr, which are nan. This is because mutations are not provided.
        - In the third and fourth column, one can see two peptides that have the same sequence but come from different vertices. Moreover, they both have the same reading frame. This is because the sequence has a stop codon in exon 1 if reading frame 0 is used, so that in both cases only the part of exon 1 up to the mutation is translated.
        - On the other hand, on the fifth line, we can see how a shift of 1 nucleotide in the reading frame leads to the disappearance of the stop codon, so that the whole sequence is translated.


Build mode *with* mutations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the second part of the example, we introduce somatic and germline mutations in the analysis. The command used in this tutorial to run the *build* mode is:

Command
~~~~~~~
.. code-block:: console

    immunopepper  build --output-dir immunopepper_usecase/ --ann-path  immunopepper/tests/data_simulated/data/build_mode/simulated_Ipp.gtf --splice-path  immunopepper/tests/data_simulated/data/build_mode/genes_graph_conf3.merge_graphs.pickle --ref-path  immunopepper/tests/data_simulated/data/build_mode/genome.fa --kmer 9 --count-path immunopepper/tests/data_simulated/data/build_mode/genes_graph_conf3.merge_graphs.count.hdf5 --parallel 1 --batch-size 1  --start-id 0 --process-num 0 --output-fasta --somatic immunopepper/tests/data_simulated/data/build_mode/variants_somatic.vcf --germline immunopepper/tests/data_simulated/data/build_mode/variants_germline.vcf --mutation-sample simulated_Ipp_1_sample3 --verbose 2

In this command, the build mode of immunopepper is run on the :ref:`input_data` described in the section above. Moreover, the output directory is set to a folder called *immunopepper_usecase/cohort_mutsimulated_Ipp_1_sample3*, located on the directory where the command is executed. The kmer length is set to 9, as it is a common kmer length selected in clinical applications. Finally, there are also two mutation files provided, a somatic and a germline file. These files will apply the existing mutations and take them into account when computing the output.
One important thing to note is that, if mutations are provided, an extra filter layer is included. This layer will ensure that only peptides different to the reference (base genome + germline) are included in the output.

Terminal output
~~~~~~~~~~~~~~~

The output displayed in the terminal is the following:

.. code-block:: console

    2023-06-20 19:21:51,580 INFO     Command lineNamespace(output_dir='immunopepper_usecase/', ann_path='immunopepper/tests/data_simulated/data/build_mode/simulated_Ipp.gtf', splice_path='immunopepper/tests/data_simulated/data/build_mode/genes_graph_conf3.merge_graphs.pickle', ref_path='immunopepper/tests/data_simulated/data/build_mode/genome.fa', kmer=9, libsize_extract=False, all_read_frames=False, count_path='immunopepper/tests/data_simulated/data/build_mode/genes_graph_conf3.merge_graphs.count.hdf5', output_samples=[], heter_code=0, compressed=True, parallel=1, batch_size=1, pickle_samples=[], process_chr=None, complexity_cap=None, genes_interest=None, start_id=0, process_num=0, skip_annotation=False, keep_tmpfiles=False, libsize_path=None, output_fasta=True, force_ref_peptides=False, filter_redundant=False, kmer_database=None, gtex_junction_path=None, disable_concat=False, disable_process_libsize=False, mutation_sample='simulated_Ipp_1_sample3', germline='immunopepper/tests/data_simulated/data/build_mode/variants_germline.vcf', somatic='immunopepper/tests/data_simulated/data/build_mode/variants_somatic.vcf', sample_name_map=None, use_mut_pickle=False, verbose=2)
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

Output files
~~~~~~~~~~~~

1. **somatic_and_germline_annot_peptides.fa.gz**: This is a fasta file containing the background peptides for each gene transcript. The file contains a header that is the transcript id and the sequence of the corresponding peptide. The name also shows the mutation mode, which in this case is somatic and germline. *Note*: As this genome is simulated, there is a higher frequency of stop codons than in nature, that explains the existence of some short peptides.

    *Output example:*

    .. code-block::

        >gene8.t1
        TSSRTMETLVP
        >gene4.t2
        SVRRTPRFRRTAEAPVSRSLIITHLGDGGWEP
        >gene15.t1
        TFLKKSLRLNSI
        ...

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
        ...


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
    | ...   |  ...                    | ...                   |         ...             | ...                    | ...                    | ...                    | ...                    | ...                    | ...                    |      ...                 | ...                   |    ...                 |   ...                  |   ...                 |  ...                   | ...                     |...                     | ...                    | ...                     |...                    |
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
        ...

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
       | ...       |  ...                            |   ...           |  ...              |       ...          |    ...               |  ...                 |   ...                |  ...                 | ...                  | ...                  | ...                  | ...                  |   ...                |   ...                 |  ...                  | ...                   |...                    | ...                   | ...                   |...                    | ...                   | ...                   | ...                   | ...                   |
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
       | ...       |  ...                            |   ...           |  ...              |       ...          |    ...               |  ...                 |   ...                |  ...                 | ...                  | ...                  | ...                  | ...                  |   ...                |   ...                 |  ...                  | ...                   |...                    | ...                   | ...                   |...                    | ...                   | ...                   | ...                   | ...                   |
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
    | ...                                                                                        | ...                           | ...                 | ...                 | ...      | ...        | ...        | ...                   | ...          | ...               | ...        | ...          | ...             | ...                 | ...                 | ...       | ...          |
    +--------------------------------------------------------------------------------------------+-------------------------------+---------------------+---------------------+----------+------------+------------+-----------------------+--------------+-------------------+------------+--------------+-----------------+---------------------+---------------------+-----------+--------------+

    Important things to note:

        - In this case, the mutation index in the id is not always 0, since we have a mutation file and we study several mutations.
        - The field *variantComb* shows the variant combination used to create the peptide. Variants in the software are SNP at DNA level.
        - The field *variantSegExpr* shows the expression of the gene with the variant combination used to create the peptide.
        - We can see once more that the sequence for the gene4 peptides is the same. This is because there is a stop codon in the exon 0 of the gene, so that the peptide is truncated.

.. _tutorial_samplespecif:

Running samplespecif mode
--------------------------

In this case, we will apply this mode to the output of build mode, in order to remove the background kmers from the foreground samples.

Command
^^^^^^^

.. code-block:: console

    immunopepper samplespecif --annot-kmer-files immunopepper_usecase/cohort_mutNone/ref_annot_kmer.gz --output-dir immunopepper_usecase/samplespecif --junction-kmer-files immunopepper_usecase/cohort_mutNone/ref_graph_kmer_JuncExpr --bg-file-path immunopepper_usecase/samplespecif/bg-file.gz --output-suffix trial

Terminal output
^^^^^^^^^^^^^^^

.. code-block:: console

    immunopepper_usecase/cohort_mutNone/ref_graph_kmer_JuncExpr --bg-file-path immunopepper_usecase/samplespecif/bg-file.gz --output-suffix trial
    2023-07-10 22:05:57,679 INFO     Command lineNamespace(annot_kmer_files=['immunopepper_usecase/cohort_mutNone/ref_annot_kmer.gz'], output_dir='immunopepper_usecase/samplespecif', junction_kmer_files=['immunopepper_usecase/cohort_mutNone/ref_graph_kmer_JuncExpr'], bg_file_path='immunopepper_usecase/samplespecif/bg-file.gz', output_suffix='trial', remove_bg=False, verbose=1)
    2023-07-10 22:05:57,679 INFO     >>>>>>>>> Remove annotation: Start
    2023-07-10 22:05:57,679 INFO     ...consider annotation file:immunopepper_usecase/cohort_mutNone/ref_annot_kmer.gz
    2023-07-10 22:05:57,709 INFO     generated unique background kmer file in immunopepper_usecase/samplespecif/bg-file.gz

    2023-07-10 22:05:57,710 INFO     ...consider foreground file:['immunopepper_usecase/cohort_mutNone/ref_graph_kmer_JuncExpr/part-b7a3198e-37eb-453b-8d92-b2dd5b855d58.gz', 'immunopepper_usecase/cohort_mutNone/ref_graph_kmer_JuncExpr/part-83795625-497b-4cf8-be61-cc363197f88d.gz', 'immunopepper_usecase/cohort_mutNone/ref_graph_kmer_JuncExpr/part-ac0ec9d0-6599-4ff6-a066-7b4538007680.gz']
    2023-07-10 22:05:57,731 INFO     output bg-removed kmer file : immunopepper_usecase/samplespecif/ref_graph_kmer_JuncExpr_trial.gz

    2023-07-10 22:05:57,731 INFO     >>>>>>>>> Remove annotation: Finish

Output files
^^^^^^^^^^^^

For the purpose of this example, we provided an empty path in `--bg-file-path` option. Therefore, this mode will also generate the intermediate file containing the unique set of background peptides. Moreover, as `--remove-bg` is set to False, the mode will generate a file identical to **ref_graph_kmer_JuncExpr** from build mode, but with an extra column *is_neo_flag* indicating whether the kmer is a neoantigen or not i.e. If the kmer was present in the background set. Finally, as we set the suffix to be *trial*, the output file will be named **ref_graph_kmer_JuncExpr_trial.gz**.

1. **ref_graph_kmer_JuncExpr_trial.gz**: This file contains the kmers from the splice graph, with an extra column *is_neo_flag* indicating whether the kmer is a neoantigen or not.

       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-------------+
       | kmer      | coord                           | isCrossJunction | junctionAnnotated | readFrameAnnotated | simulatedIpp1sample1 | simulatedIpp1sample2 | simulatedIpp1sample3 | simulatedIpp1sample4 | simulatedIpp1sample5 | simulatedIpp1sample6 | simulatedIpp1sample7 | simulatedIpp1sample8 | simulatedIpp1sample9 | simulatedIpp1sample10 | simulatedIpp1sample11 | simulatedIpp1sample12 | simulatedIpp1sample13 | simulatedIpp1sample14 | simulatedIpp1sample15 | simulatedIpp1sample16 | simulatedIpp1sample17 | simulatedIpp1sample18 | simulatedIpp1sample19 | simulatedIpp1sample20 | is_neo_flag |
       +===========+=================================+=================+===================+====================+======================+======================+======================+======================+======================+======================+======================+======================+======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=======================+=============+
       | SSSLVSDGW | 1140:1150:1300:1317:None:None   |   True          |  False            |      True          |    92.0              |  183.0               |  175.0               |  119.0               | 285.0                | 153.0                | 85.0                 | 170.0                |  40.0                |  54.0                 |  95.0                 |  99.0                 | 10.0                  |  93.0                 |  24.0                 |  87.0                 |  41.0                 |  172.0                |  69.0                 | 32.0                  | True        |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-------------+
       | EPPTYGRPSV| 1383:1400:1500:1510:None:None   |   True          |  False            |      False         |    54.0              |  117.0               |  75.0                |  77.0                | 138.0                | 70.0                 | 33.0                 | 93.0                 |  23.0                |  31.0                 |  46.0                 |  75.0                 | 6.0                   |  42.0                 |  4.0                  |  61.0                 |  26.0                 |  83.0                 |  30.0                 | 6.0                   | True        |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-------------+
       |RESSSLVSD  | 1134:1150:1300:1311:None:None   |   True          |  False            |      True          |    92.0              |  183.0               |  175.0               |  119.0               | 285.0                | 153.0                | 85.0                 | 170.0                |  40.0                |  54.0                 |  95.0                 |  99.0                 | 10.0                  |  93.0                 |  24.0                 |  87.0                 |  41.0                 |  172.0                |  69.0                 | 32.0                  | True        |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-------------+
       | ...       |  ...                            |   ...           |  ...              |       ...          |    ...               |  ...                 |   ...                |  ...                 | ...                  | ...                  | ...                  | ...                  |   ...                |   ...                 |  ...                  | ...                   |...                    | ...                   | ...                   |...                    | ...                   | ...                   | ...                   | ...                   | ...         |
       +-----------+---------------------------------+-----------------+-------------------+--------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-----------------------+-------------+

2. **bg-file.gz**: This file contains the set of unique kmers found in the background file provided.

    .. code-block:: console

        kmer
        VRSSGTSKW
        GYKFPLSSD
        RRVTNSLIM
        RHYTLGIIA
        ...

.. _tutorial_cancerspecif:

Running cancerspecif mode
--------------------------

Simulating normal data
^^^^^^^^^^^^^^^^^^^^^^

The *cancerspecif* mode takes as input the kmer files derived from the normal and cancer samples. As this example is based on simulated data, we compute the build mode on our splice graph and count information, which will correspond to the cancer data.

Therefore, we will then simulate the normal data from the cancer data. This can be done by running the code *generate_normal_data.py*. The intuition behind the simulation is as follows:

1. Drop some kmers appearing in the cancer data. These dropped kmers will be present in the cancer files but not in the normal files, the first condition for a kmer to be considered a neopeptide. Each kmer is dropped with a probability of 30%.
2. Once the kmers are dropped, we will have the full set of normal peptides. The next thing to do is to change the expression level of the kmers.

    - For some of the kmers, the expression across all the samples will be set to 0. If this is the case, the kmer will be considered as part of the annotation only and it will be removed from the final normal kmer set.
    - For the rest of the kmers, the expression is set to a random number between 0 and 105 for every sample.

3. Once the expression is set, we will have the full set of normal peptides that will be used for the filtering.

Running cancerspecif mode
^^^^^^^^^^^^^^^^^^^^^^^^^^
In this mode, the user can perform different filtering steps to keep only the kmers that are specific to a cancer sample or a cancer cohort.

For this example, the filtering will be performed in the files generated *without* mutations. Moreover, junction kmers of both cancer and normal samples will be used. There are different thresholds and filters used in this example. The different options selected will be explained below:

1. **Normal filtering:** For the normal filtering, three filters were used:

    - Filter on presence in sample: The kmers that are only present in the annotation but not in the samples are removed from the normal kmer set.
    - Filter on the number of samples: *--n-samples-lim-normal* = 15. This filter checks in how many samples (out of 20) each kmer is expressed. If it is expressed in 15 or more samples, it will be considered as a normal kmer and it will be added to the final normal kmer set.
    - Filter on the expression: *--cohort-expr-support-normal* = 100. This filter checks the expression of the kmer across all the samples. If the expression is higher than 100 in at least one sample, the kmer will be considered as a normal kmer and it will be added to the final normal kmer set.

2. **Cancer filtering:** For the cancer filtering, three types of filtering were used:

    - NeoJunctions filtering: By setting *--filterNeojuncCoord* = C, the cancer kmers are first filtered by JunctionAnnotated. Only those kmers where JunctionAnnotated = False will be considered, since we are interested in kmers resulting from novel splicing junctions.
    - Sample based filtering: We select a sample of interest *--ids-cancer-samples** = simulated_Ipp_1_sample3. This sample will be subjected to an individual filtering on its expression. The filter is set with the argument *--sample-expr-support-cancer* = 20. This filter checks the expression of the kmer in the sample of interest. If the expression is higher than 20, the kmer will be considered as a cancer kmer and it will be added as a possible cancer kmer.
    - Cohort based filtering: For the rest of the samples, excluding the sample of interest, we will perform a number of samples and expession filtering. The values are set with *--n-samples-lim-cancer* = 2 and *--cohort-expr-support-cancer* = 110. This filter checks in how many samples (out of 20) each kmer is expressed with an expression level higher than 110. Those kmers expressed in 2 or more samples with an expression level >= 110 will be considered as cancer kmers and they will be added as a possible cancer kmer.

Finally, a differential filtering between cancer kmers and normal kmers is performed. We will consider cancer kmers the ones fulfilling all the conditions (intersection of the three filters) and normal kmers the ones fulfilling at least one of the normal conditions (union of the two last filters). In the end, only the kmers belonging to the cancer set, and not present in the normal set, will be kept.

Command
~~~~~~~

The command to run the *cancerspecif* mode is:

.. code-block:: console

    immunopepper cancerspecif --cores 2 --mem-per-core 2048 --parallelism 3 --kmer 9 --output-dir immunopepper_usecase/filter_case/ --interm-dir-norm immunopepper_usecase/filter_case --interm-dir-canc immunopepper_usecase/filter_case --ids-cancer-samples simulated_Ipp_1_sample3 --mut-cancer-samples ref --output-count immunopepper_usecase/filter_case/output-count.txt --path-normal-matrix-edge immunopepper/tests/data_simulated/data/cancerspecif_mode/ref_graph_kmer_NormalExpr/normal_junctions.gz --n-samples-lim-normal 15 --cohort-expr-support-normal 100 --sample-expr-support-cancer 20 --cohort-expr-support-cancer 110 --n-samples-lim-cancer 2 --path-cancer-matrix-edge immunopepper/tests/data_simulated/data/cancerspecif_mode/ref_graph_kmer_CancerExpr/cancer_junctions.gz --filterNeojuncCoord C --verbose 1

Terminal output
~~~~~~~~~~~~~~~

The terminal output for this mode is:

.. code-block:: console

    2023-06-26 12:20:51,722 INFO     Command lineNamespace(cores=2, mem_per_core=2048, parallelism=3, out_partitions=None, scratch_dir='', interm_dir_norm='immunopepper_usecase/filter_case', interm_dir_canc='immunopepper_usecase/filter_case', kmer='9', ids_cancer_samples=['simulated_Ipp_1_sample3'], mut_cancer_samples=['ref'], whitelist_normal=None, whitelist_cancer=None, path_cancer_libsize=None, path_normal_libsize=None, normalizer_cancer_libsize=None, normalizer_normal_libsize=None, output_dir='immunopepper_usecase/filter_case/', output_count='immunopepper_usecase/filter_case/output-count.txt', tag_normals='', tag_prefix='', path_normal_matrix_segm=None, path_normal_matrix_edge=['immunopepper/tests/data_simulated/data/cancerspecif_mode/ref_graph_kmer_NormalExpr/normal_junctions.gz'], n_samples_lim_normal=15, cohort_expr_support_normal=100.0, sample_expr_support_cancer=20.0, cohort_expr_support_cancer=110.0, n_samples_lim_cancer=2, path_cancer_matrix_segm=None, path_cancer_matrix_edge=['immunopepper/tests/data_simulated/data/cancerspecif_mode/ref_graph_kmer_CancerExpr/cancer_junctions.gz'], cancer_support_union=False, path_normal_kmer_list=None, uniprot=None, filterNeojuncCoord='C', filterAnnotatedRF='', tot_batches=None, batch_id=None, on_the_fly=False, verbose=1)
    driver_mem 3072
    memory_per_executor_mb 80% 1638
    parallelism_ 3
    shuffle_partitions 3
    permsize 1024M
    2023-06-26 12:21:06,642 INFO

     >>>>>>>> Preprocessing libsizes
    2023-06-26 12:21:06,642 INFO     At least one intermediate normals filtering file is missing.
    2023-06-26 12:21:06,642 INFO     Will compute full filtering steps according to user input parameters
    2023-06-26 12:21:06,642 INFO

     >>>>>>>> Preprocessing Normal samples
    2023-06-26 12:21:06,642 INFO     Load input ['immunopepper/tests/data_simulated/data/cancerspecif_mode/ref_graph_kmer_NormalExpr/normal_junctions.gz']
    2023-06-26 12:21:12,558 INFO     ...partitions: 1
    2023-06-26 12:21:12,698 INFO     ...partitions: 1
    2023-06-26 12:21:12,698 INFO     Isolating kmers only in backbone annotation
    2023-06-26 12:21:13,043 INFO     >>>> Save to immunopepper_usecase/filter_case/kmers_derived_solely_from_annotation.tsv.gz
    2023-06-26 12:21:13,880 INFO     Cast types
    2023-06-26 12:21:14,101 INFO     ...partitions: 1
    2023-06-26 12:21:14,101 INFO

     >>>>>>>> Normals: Perform Hard Filtering
     (expressed in 15 samples with 100.0 normalized counts
    2023-06-26 12:21:14,102 INFO     expression filter
    2023-06-26 12:21:14,102 INFO     Filter matrix with cohort expression support >= 100.0 in 1 sample
    2023-06-26 12:21:14,540 INFO     Save intermediate 1/2 normals filtering file to immunopepper_usecase/filter_case/interm_normals_combiExprCohortLim100.0Across1.tsv.gz
    2023-06-26 12:21:15,655 INFO     Filter matrix with cohort expression support > 0.0 in 1 sample
    2023-06-26 12:21:16,001 INFO     Save intermediate 2/2 normals filtering file to immunopepper_usecase/filter_case/interm_normals_combiExprCohortLim0.0Across1.tsv.gz
    2023-06-26 12:21:16,615 INFO     Filter matrix with cohort expression support > 0 in 15 sample(s)
    2023-06-26 12:21:16,665 INFO     Load input immunopepper_usecase/filter_case/kmers_derived_solely_from_annotation.tsv.gz
    2023-06-26 12:21:16,976 INFO     At least one intermediate cancer_ref filtering file is missing.
    2023-06-26 12:21:16,977 INFO     Will compute full filtering steps according to user input parameters
    2023-06-26 12:21:16,977 INFO

     >>>>>>>> Preprocessing Cancer sample simulated_Ipp_1_sample3
    2023-06-26 12:21:16,977 INFO     Load input ['immunopepper/tests/data_simulated/data/cancerspecif_mode/ref_graph_kmer_CancerExpr/cancer_junctions.gz']
    2023-06-26 12:21:17,263 INFO     ...partitions: 1
    2023-06-26 12:21:17,311 INFO     ...partitions: 1
    2023-06-26 12:21:17,341 INFO     Cast types
    2023-06-26 12:21:17,473 INFO     ...partitions: 1
    2023-06-26 12:21:17,609 INFO     ...partitions: 1
    2023-06-26 12:21:17,609 INFO     Filter with simulatedIpp1sample3 > 0
    2023-06-26 12:21:18,270 INFO     # Init_cancer n = 59 kmers
    2023-06-26 12:21:18,348 INFO     ...partitions: 1
    2023-06-26 12:21:18,348 INFO     Filter with simulatedIpp1sample3 >= 20.0
    2023-06-26 12:21:18,660 INFO     # Filter_Sample n = 59 kmers
    2023-06-26 12:21:18,660 INFO     >>>> Save to immunopepper_usecase/filter_case/condition2
    2023-06-26 12:21:18,931 INFO     Target sample simulatedIpp1sample3 not included in the cohort filtering
    2023-06-26 12:21:18,931 INFO     Filter matrix with cohort expression support >= 110.0 in 1 sample
    2023-06-26 12:21:19,196 INFO     Save intermediate 1/2 cancer_ref filtering file to immunopepper_usecase/filter_case/interm_cancer_ref_combiExprCohortLim110.0Across1ExceptsimulatedIpp1sample3.tsv.gz
    2023-06-26 12:21:19,375 INFO     Filter matrix with cohort expression support > 0.0 in 1 sample
    2023-06-26 12:21:19,650 INFO     Save intermediate 2/2 cancer_ref filtering file to immunopepper_usecase/filter_case/interm_cancer_ref_combiExprCohortLim0.0Across1ExceptsimulatedIpp1sample3.tsv.gz
    2023-06-26 12:21:19,889 INFO     Filter matrix with cohort expression support >= 110.0 in 2 sample(s)
    2023-06-26 12:21:19,994 INFO     support intersect
    2023-06-26 12:21:20,778 INFO     # Filter_Sample_Cohort n = 13 kmers
    2023-06-26 12:21:20,778 INFO

     >>>>>>>> Cancers: Perform differential filtering
    2023-06-26 12:21:21,071 INFO     partitions: 1
    2023-06-26 12:21:21,071 INFO     Filtering normal background
    2023-06-26 12:21:21,822 INFO     partitions: 1
    2023-06-26 12:21:21,822 INFO     >>>> Save to immunopepper_usecase/filter_case/simulated_Ipp_1_sample3_ref_SampleLim20.0CohortLim110.0Across2_FiltNormalsCohortlim100.0Across15.tsv.gz
    2023-06-26 12:21:23,072 INFO     # Filter_Sample_Cohort_CohortNormal n = 3 kmers
    2023-06-26 12:21:24,044 INFO     Save intermediate info to immunopepper_usecase/filter_case/output-count.txt
    2023-06-26 12:21:24,879 INFO     Closing down clientserver connection

Output files
~~~~~~~~~~~~~

In this mode, there are intermediate files generated for cancer and normal samples. These intermediate files will help speed up re-runs, as filtering steps can be computationally expensive for large cohorts.

**Normal samples**

For normal samples, only the intermediate files will be generated. These files will be later used to exclude some kmers as cancer candidates, as they will not be cancer-specific kmers.

1. **kmers_derived_solely_from_annotation.tsv.gz**: This is a folder containing the kmers that are only present in the annotation and not in the samples. This means that these kmers have a zero expression in all the samples, and therefore will be excluded from the normal dataset.


    **Format:** The files contain a header, "kmer", and then a kmer in each row.

    .. code-block:: console

        kmer
        FLGVLPNAY
        LPFRVLIIL
        GVCTLGILR
        NCRKGFLGV
        GFLGVLPNA
        SSSLVSDGW
        ...

2. **interm_normals_combiExprCohortLim0.0Across1.tsv.gz:** This folder contains the kmers that are present with an expression bigger than zero in at least one sample. It is an intermediate file that will be later used for the filter based on the number of samples. The rest of the filtering, in which *--n-lim-samples-normal* threshold is applied will be performed on the fly.

    **Format:** The files obtained in this folder will be tab-separated, and they will have two columns, with the first column showing the kmer and the second column showing the number of samples in which that kmer appears with more expression than 0.

    .. code-block:: console

        KGFLGVCTL       19
        LEPATSLDF       20
        GFLGVCTLG       19
        RKGFLGVCT       20
        VLPNAYALI       19
        PFRVLIILE       20
        NCRKGFLGV       19
        KGFLGVLPN       20
        ...

3. **interm_normals_combiExprCohortLim100.0Across1.tsv.gz:** This folder contains the kmers that appear with an expression >= 100 in at least one sample. This is the file that will be directly used for the expression based filtering. The kmers passing the threshold will be marked as normal kmers and will be removed from the possible cancer kmers.

    **Format:** The file is a tab seperated file with two columns, with the first column showing the kmer and the second column showing the number of samples in which that kmer appears with more or equal expression than 100.

    .. code-block:: console

        LEPATSLDF       1
        GFLGVCTLG       3
        RKGFLGVCT       2
        VLPNAYALI       3
        PFRVLIILE       1
        NCRKGFLGV       1
        ...

**Cancer samples**

1. **interm_cancer_ref_combiExprCohortLim0.0Across1ExceptsimulatedIpp1sample3.tsv.gz:** This folder will contain the intermediate files showing the kmers having an expression bigger than 0.0 in at least 1 sample, without taking into account the target sample.

    **Format:** The files obtained in this folder will be tab separated, and they will have two columns, with the first column showing the kmer and the second column showing the number of samples in which that kmer appears with more expression than 0.

    .. code-block:: console

        KGFLGVCTL       19
        VCTLGILRV       19
        FLGVCTLGI       19
        LEPATSLDF       19
        GFLGVCTLG       19
        FLGVLPNAY       19
        ...

2. **interm_cancer_ref_combiExprCohortLim110.0Across1ExceptsimulatedIpp1sample3.tsv.gz:** This folder will contain the intermediate files showing the kmers that have an expression level higher than the expression threshold, 110 in this case, in at least one sample, without taking into account the target sample.

    **Format:** The files obtained in this folder will be tab separated, and they will have two columns, with the first column showing the kmer and the second column showing the number of samples in which that kmer appears with more expression than 110.

    .. code-block::

        SSSLVSDGW       6
        EPTYGRPSV       2
        RESSSLVSD       6
        SRHRESSSL       6
        ...

3. **simulated_Ipp_1_sample3_ref_SampleLim20.0CohortLim110.0Across2_FiltNormalsCohortlim100.0Across15.tsv.gz:** This is the file containing the cancer specific kmers, after application of the different thresholds and differential filtering against normal kmers. If an external database is not provided under *--uniprot*, this is the final output of the cancerspecif mode.

    **Format**: The output is a tab separated file made up of 5 different columns. The first column contains the cancer kmer, the second column contains the coordinates of the kmer, the third column contains the expression level of the kmer in the sample of interest, the fourth column shows whether the junction is annotated or if it is a novel junction, and the last column shows whether the read frame is annotated or not.

    +---------------+---------------------------------+----------------------+-------------------------+-------------------------+
    | kmer          |   coord                         | simulatedIpp1sample3 |   junctionAnnotated     |     readFrameAnnotated  |
    +===============+=================================+======================+=========================+=========================+
    | HRESSSLVS     |  1131:1150:1300:1308:None:None  |     175.0            |  False                  |     True                |
    +---------------+---------------------------------+----------------------+-------------------------+-------------------------+
    | ESSSLVSDG     |   1137:1150:1300:1314:None:None |     175.0            |  False                  |     True                |
    +---------------+---------------------------------+----------------------+-------------------------+-------------------------+
    | SLVSDGWAC     |  1146:1150:1300:1323:None:None  |     175.0            |  False                  |     True                |
    +---------------+---------------------------------+----------------------+-------------------------+-------------------------+



4. **output-count.txt:** This file will be generated because *--output-count* was provided in the arguments. It contains the number of remaining kmers after each filtering step for each sample of interest.

    +------------------------+-----------------+--------------------+--------------------------+----------------------------+-----------------------------+----------------------+--------------------------------+---------------+-------------------+-----------------------+------------------------------------+
    |sample                  |    mutation_mode|   min_sample_reads |       #_of_cohort_samples|     reads_per_cohort_sample|   #_normal_samples_allowed  |      normal_cohort_id|        reads_per_normal_sample |   Init_cancer |    Filter_Sample  | Filter_Sample_Cohort  |   Filter_Sample_Cohort_CohortNormal|
    +========================+=================+====================+==========================+============================+=============================+======================+================================+===============+===================+=======================+====================================+
    |simulated_Ipp_1_sample3 |    ref          |   20.0             |       2                  |     110.0                  |   15                        |                      |        100                     |   59          |    59             | 13                    |   3                                |
    +------------------------+-----------------+--------------------+--------------------------+----------------------------+-----------------------------+----------------------+--------------------------------+---------------+-------------------+-----------------------+------------------------------------+


.. _tutorial_mhcbind:

Running mhcbind mode
---------------------

This mode is a wrapper tool for the mhctools package. It allows to do mhc binding predictions on the output kmers from *cancerspecif* mode. Therefore, there are some input parameters specific to this software and some input parameters that will depend on the mhc prediction tool that the user chooses to use.

In this example, a binding affinity prediction on the three kmers given as output from the *cancerspecif* mode will be done. The affinity prediction will be done using the tool *mhcflurry*, and a threshold for affinity >=1000 is set to select the final kmers.

.. note:: The mhctools package requires the mhc binding affinity tools to be locally installed in the user's computer. For more information on how to install the mhc binding affinity tools, please refer to the mhctools documentation.

Command
^^^^^^^^

The command to run the *mhcbind* mode is the following:

.. note:: In this case you need to set the --mhc-software-path argument to the path where your cloned "mhctools" repository is.

.. code-block::

    immunopepper mhcbind --mhc-software-path ./immunopepper/mhctools/ --argstring "--mhc-predictor mhcflurry --mhc-alleles HLA-A*02:01 --output-csv immunopepper_usecase/mhc_bind/predictions.csv --input-peptides-file immunopepper_usecase/mhc_bind/input_peptides.csv" --partitioned-tsv immunopepper_usecase/filter_case/simulated_Ipp_1_sample3_ref_SampleLim20.0CohortLim110.0Across2_FiltNormalsCohortlim100.0Across15.tsv.gz --output-dir immunopepper_usecase/mhc_bind --bind-score-method affinity --bind-score-threshold 1000 --verbose 1

All the arguments *outside* the --argstring belong to the immunopepper software, while the ones inside the argstring are used by the mhctools.

**Immunopepper arguments:**

- *--mhc-software-path*: Path to the mhctools package. This is a required argument.
- *--partitioned-tsv*: Path to the output folder of the *cancerspecif* mode. If the user wants to run the *mhcbind* mode on the output of the *cancerspecif* mode, this argument is required.
- *output-dir*: Path to the output folder. This is a required argument. It should match the directory given under *--output-csv* in the *--argstring* argument.
- *--bind-score-method*: Metric used for filtering. The filtering will be done on the results of the mhc binding prediction, so the score method should be an output of the selected mhc binding prediction tool.
- *--bind-score-threshold*: Threshold for the filtering. The filtering will be done on the results of the selected metric in the previous argument. The threshold should be a number.
- *--verbose*: Verbosity level. This is an optional argument. The default value is 1.

**Mhctools arguments:**

The arguments for mhctools are the ones contained in the *--argstring*:

- *--mhc-predictor*: Name of the mhc binding prediction tool to use.
- *--output-csv*: Path to the output file and name of the output file. The format should be *.csv*.
- *--input-peptides-file*: Path to the input file containing the kmers to predict. The kmers are derived from the file provided under *--partitioned-tsv*. The format should be *.csv*.
- *--mhc-alleles*: Alleles to use for the prediction.

Terminal output
^^^^^^^^^^^^^^^

The output displayed in the terminal for the *mhcbind* mode is the following:

.. code-block:: console

    2023-06-27 15:53:05,501 INFO     Command lineNamespace(mhc_software_path='./immunopepper/mhctools/', argstring='--mhc-predictor mhcflurry --mhc-alleles HLA-A*02:01 --output-csv immunopepper_usecase/mhc_bind/predictions.csv --input-peptides-file immunopepper_usecase/mhc_bind/input_peptides.csv', output_dir='immunopepper_usecase/mhc_bind', partitioned_tsv='immunopepper_usecase/filter_case/simulated_Ipp_1_sample3_ref_SampleLim20.0CohortLim110.0Across2_FiltNormalsCohortlim100.0Across15.tsv.gz', bind_score_method='affinity', bind_score_threshold=1000.0, less_than=False, verbose=2)
    2023-06-27 15:53:05,892 INFO     Process the outputs from cancerspecif mode
    2023-06-27 15:53:05,905 INFO     Launch MHC Tools with command --mhc-predictor mhcflurry --mhc-alleles HLA-A*02:01 --output-csv immunopepper_usecase/mhc_bind/predictions.csv --input-peptides-file immunopepper_usecase/mhc_bind/input_peptides.csv
    2023-06-27 15:53:05,906 - mhctools.cli.args - INFO - Building MHC binding prediction type for alleles ['HLA-A*02:01'] and epitope lengths None
    2023-06-27 15:53:05,906 INFO     Building MHC binding prediction type for alleles ['HLA-A*02:01'] and epitope lengths None
    2023-06-27 15:53:11,257 INFO     Loaded 10 class1 pan allele predictors, 14839 allele sequences, 6308 percent rank distributions, and 0 allele specific models:
    Wrote: immunopepper_usecase/mhc_bind/predictions.csv
    2023-06-27 15:53:21,936 INFO     Perform filtering with affinity >= 1000.0
    2023-06-27 15:53:21,938 INFO     Saving to immunopepper_usecase/mhc_bind/predictions_WithaffinityMoreLim1000.0.tsv


Output files
^^^^^^^^^^^^

This mode creates files that are unique to *immunopepper* or files that are created by the mhctools package. This will be indicated in the output file description.

1. **input_peptides.csv**: This is the file created from *--partitioned-tsv*. The name is given by the user, when setting the argument *--input-peptides-file* in the argstring. It is a file unique to *immunopepper*, as it is the intermediate file that makes compatible the output of *cancerspecif* mode and the input to mhctools.

    **Format**: *.csv* file with one column and without header. It contains a kmer per row. The kmers are the ones that passed the filtering steps in the *cancerspecif* mode.

    .. code-block::

        HRESSSLVS
        ESSSLVSDG
        SLVSDGWAC

2. **predictions.csv**: This is the output of the *mhc binding tool* selected. In this case, the selected tool was *mhcflurry*. This is the file generated by the *mhctools* package. The name is selected by the user when setting *--output-csv* in the argstring.

    **Format**: The format will depend on the prediction tool selected. The different prediction scores are the ones that can be used in *--bind-score-method*. In this case the options would be: score, affinity, percentile_rank.

    +---------------------+--------+-------------+--------------+-----------------+--------------------+-----------------------+-------------------------------+----------+
    |source_sequence_name |  offset|   peptide   |       allele |   score         |         affinity   |       percentile_rank |        prediction_method_name |   length |
    +=====================+========+=============+==============+=================+====================+=======================+===============================+==========+
    | None                |       0| HRESSSLVS   | HLA-A*02:01  | 0.0446471281349 | 30844.2694544      | 49.96275              | mhcflurry                     |        9 |
    +---------------------+--------+-------------+--------------+-----------------+--------------------+-----------------------+-------------------------------+----------+
    |None                 |       0| ESSSLVSDG   | HLA-A*02:01  | 0.0493709236383 | 29307.415067       | 36.933125             | mhcflurry                     |        9 |
    +---------------------+--------+-------------+--------------+-----------------+--------------------+-----------------------+-------------------------------+----------+
    |None                 |       0| SLVSDGWAC   | HLA-A*02:01  | 0.367525291443  | 937.518149464      | 1.700125              | mhcflurry                     |        9 |
    +---------------------+--------+-------------+--------------+-----------------+--------------------+-----------------------+-------------------------------+----------+

3. **predictions_WithaffinityMoreLim1000.0.tsv:** This is the output of the *mhcbind* mode. This output combines the mhc binding predictions of *mhctools* with the output of *cancerspecif* mode. Moreover, it performs filteirng if selected. In this case, filtering was selected, so the name of the file reflects the metric and the threshold used. By looking at *predictions.csv* one can see that the affinity of the last kmer is 937.518149464. This is below the threshold of 1000, so this kmer is not included in the output file.

    **Format**: The file will be in *.tsv* format. It will contain the columns of *predictions.csv*, merged with the output of *cancerspecif* mode.


   +------------+------------------------+----------------------+------------------------+--------------------------+----------+------------+---------------+---------------+---------------------+-------------------------+--------+
   | kmer       |   simulatedIpp1sample3 |   junctionAnnotated  |     readFrameAnnotated |     source_sequence_name |   offset | allele     |  score        |    affinity   |     percentile_rank |   prediction_method_name|  length|
   +============+========================+======================+========================+==========================+==========+============+===============+===============+=====================+=========================+========+
   | HRESSSLVS  |                    175 | False                | True                   | None                     |        0 | HLA-A*02:01| 0.0446471     | 30844.3       | 49.9628             | mhcflurry               |      9 |
   +------------+------------------------+----------------------+------------------------+--------------------------+----------+------------+---------------+---------------+---------------------+-------------------------+--------+
   | ESSSLVSDG  |                    175 | False                | True                   | None                     |        0 | HLA-A*02:01| 0.0493709     | 29307.4       | 36.9331             | mhcflurry               |      9 |
   +------------+------------------------+----------------------+------------------------+--------------------------+----------+------------+---------------+---------------+---------------------+-------------------------+--------+





