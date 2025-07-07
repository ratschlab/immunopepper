Output files
=============

In this section, the output files of the different modes are described.

- :ref:`build_out`: Output files generated in the *build* mode.
- :ref:`samplespecif_out`: Output files generated in the *samplespecific* mode.
- :ref:`cancerspecif_out`: Output files generated in the *cancerspecific* mode.
- :ref:`mhcbind_out`: Output files generated in the *mhcbind* mode.
- :ref:`pepquery_out`: Output files generated in the *pepquery* mode.

.. _build_out:

Outputs mode **build**
-----------------------
There are 10 output files for the *build* mode.

In the output files name, **mut_mode** refers to *ref*, *somatic*, *germline* and *somatic_and_germline*. The modes are obtained by passing the respective --germline and --somatic files to the command line. The default behaviour without variants or mutations is the reference "ref" mode.

.. _output-1-build:

1. **\[mut_mode\]_annot_peptides.fa**: FASTA file containing background peptides for each gene transcript. Generated only if --skip-annotation is set to False.

    **Format:** The header shows the transcript ID and the value is the result peptide.

    .. code-block::

        >ENST00000514520.5
        MKLSMKNNIINTQQSFVTMPNVIVPDIEKEIRRMENGACSSFSEDDDSASTSEESENENPHARGSFSYKSLRKGGPSQREQYLPGAIALFNVNNSSNKDQ


2. **\[mut_mode\]_annot_kmers**: kmers generated from **\[mut_mode\]_annot_peptides.fa**. Generated only if --skip-annotation is set to False.

    **Format:** One column called *kmer* containing each resulting kmer.

    .. code-block::

        kmer
        WLILDYVSD
        VIIIHWNAC
        TEYPDAKTM

3. **gene_expression_detail**: File containing expression for each coding gene and each sample. Generated only if --count-path file is provided..

    .. code-block::

        _____________________________________________________
         gene                  | sample 1 | ... | sample n
        _____________________________________________________
         ENSMUSG00000025902.13 | 366155   | ... | 0
                               |          |     |
         ENSMUSG00000025903.14 | 0        | ... | 0
        _____________________________________________________


4. **\[mut_mode\]_sample_peptides.fa**: FASTA file containing foreground peptides for each gene transcript obtained from traversing the splicegraph. Generated only if --output-fasta is set to True.

    **Format:** The header shows the transcript ID and the value is the result peptide.

    .. code-block::

        >ENST00000513178.1
        MKLSMKNNIINTQQSFVTMPNVIVPDIEKEIRRMENGACSSFSEDDDSASTSEESENENPHARGSFSYKSLRKGGPSQREQYLPGAIALFNVNNSSNKD

.. _output-5-build:

5. **\[mut_mode\]_graph_kmer_SegmExpr**: Folder containing a file with the expression levels of kmers found in an exon. Kmers shown in this file are generated from a single exon and are not located across an exon junction.

    **Format:** [kmer, coord, isCrossJunction, junctionAnnotated, readFrameAnnotated, sample names (nº columns = nº samples)]

    * *kmer*: str. The kmer sequence
    * *coord*: str. The genomic coordinates of the kmer. **Format:** start:end:nan:nan:None:None
    * *isCrossJunction*: Boolean. True if the kmer is located in an exon junction. In this case, all the kmers will get False.
    * *junctionAnnotated*: Boolean. True if the junction kmer is appearing in the annotation file provided under `--ann-path`. In this case, all kmers will get False.
    * *readFrameAnnotated*: Boolean. True if the reading frame was appearing in the annotation file provided under `--ann-path`. The flag value will be False if the reading frame was obtained with the immunopepper tool.

    .. code-block::

        _________________________________________________________________________________________________________________________________________________________________
        kmer       |  coord                                  |  isCrossJunction   |  junctionAnnotated  |    readFrameAnnotated   |   sample 1  |     ...  |   sample n
        _________________________________________________________________________________________________________________________________________________________________
        SDFQLIFVF  |    47970995:47971022:nan:nan:None:None  |      False         |        False        |           True          |       0.0   |     ...  |      0.0

        TLIPEELEP  |    47948807:47948834:nan:nan:None:None  |      False         |        False        |           False         |       0.0   |     ...  |      0.25
        _________________________________________________________________________________________________________________________________________________________________

.. _output-6-build:

6. **\[mut_mode\]_graph_kmer_JuncExpr**: Folder containing a file with the expression levels of kmers located across exon junctions.

    **Format:** [kmer, coord, isCrossJunction, junctionAnnotated, readFrameAnnotated]

    * *kmer*: str. The kmer sequence
    * *coord*: str. The sorted exon coordinates of the kmer. **Format:** left_e1:right_e1:left_e2:right_e2:left_e3:right_e3, with left < right. For *+* strand this corresponds to: start_e1:end_e1_start_e2:end_e2:start_e3:end_e3. For *-* strand: end_e1:start_e1:end_e2:start_e2:end_e3:start_e3. If the kmer is only crossing one junction, start_e3 and end_e3 will be `None`.
    * *isCrossJunction*: Boolean. True if the kmer is located across an exon junction. In this case, all the kmers will get True.
    * *junctionAnnotated*: Boolean. True if the junction kmer is appearing in the annotation file provided under --ann-path. If this flag is False, it means that it is a novel kmer.
    * *readFrameAnnotated*: Boolean. True if the reading frame was appearing in the annotation file provided under --ann-path. The flag value will the False if the reading frame was obtained with the immunopepper tool.


    .. code-block::

        ___________________________________________________________________________________________________________________________________________________________________________
        kmer       |  coord                                            |  isCrossJunction   |  junctionAnnotated  |    readFrameAnnotated   |   sample 1  |     ...  |   sample n
        ___________________________________________________________________________________________________________________________________________________________________________
        KKEKKSQMI  |    47943370:47943387:47943274:47943284:None:None  |      True          |        False        |           False         |       0.0   |     ...  |      1.0

        REPEEKKKK  |    47952582:47952584:47943387:47943412:None:None  |      True          |        False        |           True          |       0.0   |     ...  |      0.25
        ___________________________________________________________________________________________________________________________________________________________________________


7. **expression_counts.libsize.tsv**: File containing 75% of expression and total expression for each sample. For a given sample, the "75% of expression" is defined as the 75th quantile of the gene expression distribution across coding genes. For a given sample the "total expression" is defined as the total gene expression across coding genes. Generated only if `--disable-libsize` is set to False and if `--count-path` file is provided.

    **Format:** sample_id, 75% expression, total expression.

    .. code-block::

        _______________________________________________________________
        sample     |   libsize_75percent    |     libsize_total_count
        _______________________________________________________________
        sample 1   |       309287.0         |          6256400944.0

        sample 2   |       197045.0         |          4167429408.0
        _______________________________________________________________


8. **Annot_IS_SUCCESS**: Empty file indicating that background generation was successful. Generated only if `--skip-annotation` is set to False.

9. **Output_sample_IS_SUCCESS**: Empty file created if the foreground generation was successful.

.. _output-10-build:

10. **\[mut_mode\]_sample_peptides_meta**: File containing details for each peptide generated from an exon pair.

    Detailed explanation for columns in **\[mut_mode\]_sample_peptides_meta**:

    - **peptide**: str. The peptide sequence for a specific exon pair of the foreground data.
    - **id**: In the format of \[gene_name\]:\[first vertex\]_\[second vertex\]:\[variant_id\]:\[translation_start_position\]:\[kmerType\]. Eg: *ENSG00000198515.13:18_15:0:47952701:2-exons*. *ENSG00000198515.13* is the gene name, *18_15* means this junction consists of vertex 18 and vertex 15. *0* means there is no somatic mutation or that it is the first case of all somatic mutation combinations. *47952701* is the translation start position, and *2-exons* means that the peptide is made by the combination of 2 exons.
    - **readFrame**: int (0,1,2). Reading frame used for translating the peptide. **Note**: The reading frame number is set internally by the software and does not map to the reading frame code from the annotation file.
    - **readFrameAnnotated**: True if the reading frame was appearing in the annotation file provided under --ann-path. The flag value will be False if the reading frame was obtained with the immunopepper tool.
    - **geneName**: str. The name of gene.
    - **geneChr**: str. The chromosome where the gene is located.
    - **geneStrand**: str (+, -). The strand of gene.
    - **mutationMode**: str (ref, somatic, germline, somatic_and_germline). Mutation mode
    - **hasStopCodon**: int. Indicates if there is a stop codon in the junction pair.
    - **isJunctionList**: (np.nan, 1, 0). Indicates if the junction pair appears in the given junction whitelist provided under --gtex-junction-path.
    - **isIsolated**: int. Indicate if the output peptide is translated from a single exon instead of two.
    - **variantComb**: If mutation files are provided, it shows the somatic mutation combination used in this line of output. eg. 5;25 means the somatic mutation of position 5 and 25 take effect in this output.
    - **variantSegExpr**: If mutation file and count file are provided, this field shows the expression of segments where the somatic mutation is in.  If mutation files are not provided the value is set to nan. eg. 257.0;123.2 means the segment where the somatic mutation is in position 5 has counts 257.0, and the segment where the somatic mutation is in position 25 has counts 123.2.
    - **modifiedExonsCoord**: Sorted coordinates for the exons forming the junction. They are obtained taking into account the CDS reading frame. Format: Usually we have 4 numbers: left_v1;right_v1;left_v2;right_v2, with left < right. For *+* strand: start_v1;end_v1;start_v2;end_v2. For *-* strand: end_v1;start_v1;end_v2;start_v2.
    - **originalExonsCoord**: Shows the original exon coordinates obtained from splice graph without taking into account the CDS.
    - **vertexIdx**: Shows the vertex id of the given junction. **Note**: Vertex numbering is relative to the splicegraph. Eg: 5,6 means this junction pair consists of the fifth and sixth vertex.
    - **kmerType**: str. Shows whether the peptide was translated from 2 or 3 exons.

    .. code-block::

        _______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
         peptide                                 |  id                                          |    readFrame    |   readFrameAnnotated   |       geneName        |     geneChr   |   geneStrand  |  mutationMode |  hasStopCodon  |  isJunctionList |   isIsolated  |    variantComb |    variantSegExpr     |       modifiedExonsCoord                  |            originalExonsCoord              |     vertexIdx  |  kmerType
        _______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
         MKLSMKNNIINTQQSFVTMPNVIVPDIEKEIRRMENGA  | ENSG00000198515.13:18_15:0:47952701:2-exons  |       2         |          True          |   ENSG00000198515.13  |      chr4     |       -       |      ref      |       0        |       nan       |       0       |       nan      |         nan           |   47952582;47952701;47951354;47951469     |     47952582;47952807;47951352;47951469    |        18;15   |   2-exons

         KSDDKNENKNDPEKKKKKKDKEKKKKEEKSKDKKEEEK  | ENSG00000198515.13:7_4:0:47943287:2-exons    |       2         |          False         |   ENSG00000198515.13  |      chr8     |       +       |      ref      |       0        |       nan       |       0       |       nan      |         nan           |   47943180;47943287;47942042;47942148     |      47943180;47943288;47942040;47942148   |         7;4    |   2-exons
        _______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


.. _samplespecif_out:

Outputs mode `samplespecif`
---------------------------

This mode can have either one or two outputs. If a file containing unique background kmers is provided under `--bg-file-path`, the mode will just return the main output file. On the other hand, if this file is not provided it will be generated and returned as a mode output.

1. **\[mut_mode\]_graph_kmer_JuncExpr/SegmExpr_{output-suffix}.gz**: This is a modified version of the files in **\[mut_mode\]_graph_kmer_JuncExpr** or **\[mut_mode\]_graph_kmer_SegmExpr**. If `--remove-bg` is set to False, the files will contain a new column called `is_neo_flag`. This flag will be True if the kmer is unique to the foreground data and False if it is also present in the background data. If `--remove-bg` is set to True, the mode will return the files without the kmers that are common with the background data.

    .. code-block::

            ___________________________________________________________________________________________________________________________________________________________________________________________
            kmer       |  coord                                            |  isCrossJunction   |  junctionAnnotated  |    readFrameAnnotated   |   sample 1  |     ...  |   sample n   | is_neo_flag |
            ___________________________________________________________________________________________________________________________________________________________________________________________
            KKEKKSQMI  |    47943370:47943387:47943274:47943284:None:None  |      True          |        False        |           False         |       0.0   |     ...  |      1.0     |      True   |

            REPEEKKKK  |    47952582:47952584:47943387:47943412:None:None  |      True          |        False        |           True          |       0.0   |     ...  |      0.25    |      False  |
            ___________________________________________________________________________________________________________________________________________________________________________________________

2. **{bg-file-path}.gz**: This file will be generated only if a non-existent file is provided under `--bg-file-path`. It will contain the unique set of background kmers from all the kmers present in `--annot-kmer-files`. It will contain a kmer per row, with *kmer* as header.

    .. code-block::

        kmer
        WLILDYVSD
        VIIIHWNAC
        TEYPDAKTM

.. _cancerspecif_out:


Outputs mode `cancerspecif`
------------------------------

In this mode, if the argument --on-the-fly is set to False, there will be intermediate files generated for cancer samples and for normal samples. These intermediate files will help speed up re-runs, as filtering steps can be computationally expensive for large cohorts.

**Output files for normal samples**

For normal samples, only the intermediate files will be generated. These files will be later used to exclude some kmers as cancer candidates, as they will not be cancer-specific kmers. These files will be stored in the scratch directory if --scratch-dir is provided, in the directory specified by the user under --interm-dir-norm or in the output directory specified under --output-dir if none of the previous arguments are not provided.

In the filtering pipeline, kmers are filtered based on their expression level and on the number of samples in which they are expressed. However, the building of the splicegraph includes kmers never seen in the RNA-samples and solely present in the annotation. These kmers will be processed separately and will be removed from the foreground set. The kmers will be stored in the file **kmers_derived_solely_from_annotation.csv**.

1. **kmers_derived_solely_from_annotation.csv**: This is a folder containing the kmers that are derived from the annotation and not present in any normal sample. The kmers present in this file will be skipped, and they will not be further processed. They will also be removed from the foreground matrix.

    **Technical note:** If the expression data is taken from junctions, provided under --path-cancer-matrix-edge, the kmers selected in this file will be the ones with *JunctionAnnotated = True* and zero expression across all normal samples. If the expression data is taken from segments, provided under --path-cancer-matrix-segm, the kmers selected in this file will be the ones with *ReadFrameAnnotated = True* and zero expression across all normal samples.

2. **interm_normal_combiExprCohortLim{--cohort-expr-support-normal}Across1_batch{--batch_id}_{--tot-batches}.tsv.gz**: This folder contains intermediate calculations for the filtering of normals. It will show the kmers present with an expression bigger than --cohort-expr-support-normal in at least one normal sample.

    **Technical note:** As the filtering data is usually large, it is divided in different parts or batches and analyzed separately. In this directory one can obtain the intermediate information for each part.

    **Format:** The file is a tab seperated file with two columns, with the first column showing the kmer and the second column showing the number of samples in which that kmer appears with more expression than the expression threshold.

    **Pipeline relevance:** The expression filter --cohort-expr-support-normal is defined as the "expression threshold" needed in at least one background sample for exclusion. Therefore, the kmers present in these files are the ones passing the expression filter threshold for background samples and will be removed from the foreground. This file will be helpful to perform the expression filtering step 2a) of :ref:`normal filtering <filt-normal>`.

    The directory also contains an empty file '_SUCCESS' that indicates that the filtering was successful.

    .. code-block::

        AAAAAAAGD      7868
        AAAAAAAKN      7876
        AAAAAAAKP       2


3. **interm_normal_combiExprCohortLim0.0Across1_batch{--batch_id}_{--tot-batches}.tsv.gz**: This folder contains intermediate calculations for the filtering of normals. The intermediate files will show the kmers present with an expression bigger than 0 in at least one normal sample.

    **Technical note:** As the filtering data is usually large, it is divided in different parts or batches and analyzed separately. In this directory one can obtain the intermediate information for each part.

    **Format:** The files obtained in this folder will be tab-separated, and they will have two columns, with the first column showing the kmer and the second column showing the number of samples in which that kmer appears with more expression than 0.

    **Pipeline relevance:** The filter --n-samples-lim-normal is defined as the "recurrence threshold with any read (>0)" needed in --n-samples-lim-normal background samples for exclusion. This file will be helpful to perform the filtering step 2b) in :ref:`normal filtering <filt-normal>`.

    **Technical note:** In the filtering step 2b), the recurrence threshold defined in **pipeline relevance** will be performed on the fly using this intermediate file as input. The operation can be performed on the fly because it requires less computational power.

    The directory also contains an empty file '_SUCCESS' that indicates that the filtering was successful.

    .. code-block::

        AAAAAAAGD      7868
        AAAAAAAKN      7876
        AAAAAAAKP       50


**Output files for cancer samples**

For cancer samples, there will be both intermediate and output files generated. The intermediate files will be later used for selecting kmer candidates based on their expression and recurrence. These files will be stored in the scratch directory if `--scratch-dir` is provided, in the directory specified by the user under `--interm-dir-cancer` or in the output directory specified under `--output-dir` if the previous arguments are not provided.

1. **Intermediate files**

   a. **interm_cancer_{mutation_mode}_combiExprCohortLim{--cohort-expr-support-cancer}Across1Except{--ids-cancer-sample}_batch{--batch_id}_{--tot-batches}.tsv.gz**: This folder contains intermediate calculations for the filtering of cancer samples. It will show the kmers present with an expression higher than --cohort-expr-support-cancer in at least one cancer sample other than the target. The target sample is provided under --ids-cancer-samples and will be referenced in the output folder name.

        **Technical note:** As the filtering data is usually large, it is divided in different parts or batches and analyzed separately. In this directory one can obtain the intermediate information for each part.

        **Format:** The file is a tab-seperated file with two columns, with the first column showing the kmer and the second column showing the number of samples in which that kmer appears with more expression than the expression threshold.

        **Pipeline relevance:** The filter --cohort-expr-support-cancer is defined as the "expression threshold" requested in the cancer cohort samples, excluding the target sample. The kmers present in this intermediate file will be the ones passing --cohort-expr-support-cancer in one or more samples. However, this expression filter needs to be combined with the number of samples filter --n-samples-lim-normal, so this file will not directly show the kmers that will be selected as cancer candidates. This file will be helpful to perform the expression filtering of step of :ref:`cancer filtering <filt-cancer>`.

        The directory also contains an empty file '_SUCCESS' that indicates that the filtering was successful.

        One example of how the output name for this file could look is the following: *interm_cancer_somatic_combiExprCohortLim3.0Across1ExceptTCGA-OR-A5J1-01A-11R-A29R-07_batch0_1.tsv.gz*.

        .. code-block::

            BBBBBBBGD      6898
            BBBBBBBKN      7356
            BBBBBBBKP       900


  b. **interm_cancer_{mutation_mode}_combiExprCohortLim0.0Across1Except{--ids-cancer-sample}_batch{--batch_id}_{--tot-batches}.tsv.gz**: This folder contains intermediate calculations for the filtering of cancer samples. It will show the kmers present with an expression bigger than 0 in at least one cancer sample other than the target. The target sample is provided under --ids-cancer-samples, and will be referenced in the output folder name.

       **Technical note:** As the filtering data is usually large, it is divided in different parts or batches and analyzed separately. In this directory one can obtain the intermediate information for each part.

       **Format:** The files obtained in this folder will be tab separated, and they will have two columns, with the first column showing the kmer and the second column showing the number of samples in which that kmer appears with more expression than 0.

       **Pipeline relevance:** The filter --n-samples-lim-normal is defined as the "recurrence threshold with --cohort-expr-support-cancer reads" requested in --n-samples-lim-normal cancer cohort samples, excluding the target sample. However, this intermediate file only assesses kmer expression being bigger than 0 in one or more samples, and the --n-samples-lim-normal threshold will be applied in step 2 of :ref:`cancer filtering <filt-cancer>`.

       **Technical note:** The filtering step 2 applies the recurrence threshold defined under pipeline relevance. It is performed on the fly using this intermediate file as input. The operation can be performed on the fly because it requires less computational power.

       The directory also contains an empty file '_SUCCESS' that indicates that the filtering was successful.

       One example of how the output name for this file could look is the following: *interm_cancer_somatic_combiExprCohortLim0.0Across1ExceptTCGA-OR-A5J1-01A-11R-A29R-07_batch0_1.tsv.gz*.

       .. code-block::

            BBBBBBBGD      3056
            BBBBBBBKN      2576
            BBBBBBBKP       900



2. **Output files**

    .. _output-tsv-cancerspecif:

    a. **{tag_prefix}_{id_cancer_sample}_{mutation_mode}_SampleLim{--sample-expr-support-cancer}_CohortLim{--cohort-expr-support-cancer}_Across{--n-samples-lim-cancer}_FiltNormals{--tag_normals}_CohortLim{--cohort-expr-support-normal}_Across{--n-samples-lim-normal}_batch{--batch-id}_{--tot-batches}.tsv**: This file is obtained after the "cancer/foreground support filtering" and "normal/background differential filtering" steps. The steps aim at assessing the cancer support condition, and at retaining cancer-specific kmers only. Therefore, the result contains the kmers that:

       - Passed the cancer filter based on expression in the target sample
       - Passed the cancer cohort filter based on the expression and minimum number of samples
       - Are not present in the normal files, subject to expression and sample limits set in the normals.

       **Technical note:** `batch{--batch-id}_{--tot-batches}` is added to the folder name only if the `--batch` arguments are provided.

       One example of how the output name for this file could look is the following: *breast_TCGA-OR-A5J1-01A-11R-A29R-07_somatic_SampleLim3.0_CohortLim0.0Across10_FiltNormalsGtexcoreCohort_ExceptTCGA-OR-A5J1-01A-11R-A29R-07_CohortLim0.0_Across1_batch0_1.tsv.gz*.

       .. code-block::

           kmer            coord                                  TCGAA2A0SX01A12RA08407all       JunctionAnnotated   readFrameAnnotated
          AAGDDENHN        970:990:1000:1008:None:None                   244.0                        True                 True
          AAPGQHLQA        1300:1304:1500:1567:None:None                 38.0                         False                False
          ACSNFIFKH        2045:2050:2078:2090:None:None                 114.0                        False                True

    b. **{tag_prefix}_{id_cancer_sample}_{mutation_mode}_SampleLim{--sample-expr-support-cancer}_CohortLim{--cohort-expr-support-cancer}_Across{--n-samples-lim-cancer}_FiltNormals{--tag_normals}_CohortLim{--cohort-expr-support-normal}_Across{--n-samples-lim-normal}_FiltUniprot_batch{--batch-id}_{--tot-batches}.tsv**: This file will only be generated if the user provides an Uniprot database, under `--uniprot`, containing kmers that should not appear in the output of the foreground kmers.

       **Pipeline relevance:** This file is obtained by filtering the output file a with the kmers that are present in the provided uniprot database.

       One example of how the output name for this file could look is the following: *breast_TCGA-OR-A5J1-01A-11R-A29R-07_somatic_SampleLim3.0_CohortLim0.0Across10_FiltNormalsGtexcoreCohort_ExceptTCGA-OR-A5J1-01A-11R-A29R-07_CohortLim0.0_Across1_FiltUniprot_batch0_1.tsv.gz*.


       .. code-block::

          kmer            coord                                  TCGAA2A0SX01A12RA08407all       JunctionAnnotated   readFrameAnnotated
          AAGDDENHN        970:990:1000:1008:None:None                   244.0                        True                 True
          AAPGQHLQA        1300:1304:1500:1567:None:None                 38.0                         False                False
          ACSNFIFKH        2045:2050:2078:2090:None:None                 114.0                        False                True


    .. _output-count-cancerspecif:

    c. **{--output-count}**: If --output-count is provided, this file will be generated. It contains the number of remaining kmers after each filtering step. It is a tabular file with different fields:

       - Sample: The cancer sample id. It is an id matching *--ids-cancer-samples*.
       - Mutation mode: The mutation mode used for the analysis. It is either *ref*, *somatic*, *germline* or *somatic_and_germline*.
       - Min_sample_reads: Expression threshold for sample filtering (step 2a of :ref:`cancer pipeline <filt-cancer>`). It is the threshold provided under *--sample-expr-support-cancer*.
       - Number (#) of cohort samples: Number of samples threshold for cancer cohort filtering (step 2b of :ref:`cancer pipeline <filt-cancer>`). It is the threshold provided under *--n-samples-lim-cancer*.
       - Reads per cohort sample: Expression threshold for expression cancer cohort filtering (step 2b of :ref:`cancer pipeline <filt-cancer>`). It is the threshold provided under *--cohort-expr-support-cancer*.
       - Number (#) of normal samples allowed: Number of samples threshold for normal cohort filtering (step 2b of :ref:`normal pipeline <filt-normal>`). It is the threshold provided under *--n-samples-lim-normal*.
       - Normal cohort id: The normal cohort id.
       - Reads per normal sample: Expression threshold for expression normal cohort filtering (step 2a of :ref:`normal pipeline <filt-normal>`). It is the threshold provided under *--cohort-expr-support-normal*.
       - Init_cancer: Number of kmers before any filtering step.
       - Filter_sample: Number of kmers remaining after sample filtering (step 2a of :ref:`cancer pipeline <filt-cancer>`).
       - Filter_sample_cohort: Number of kmers remaining after cancer cohort filtering (step 2b of :ref:`cancer pipeline <filt-cancer>`).
       - Filter_sample_cohort_cohortnormal Number of kmers after differential filtering with normal samples (i.e. After removing from the foreground the kmers that passed the thresholds for expression *--cohort-expr-support-normal* and number of samples *--n-samples-lim-normal* in the normal cohort).
       - Filter_sample_cohort_cohortnormal_uniprot: Number of kmers after filtering with uniprot selected database.
       - Info: If a *--tag-normals* is provided, it will be indicated here. **Technical note**: *tag-normals* can be any information specific to the normal cohort that the user wishes to record (e.g. name of the cohort or the total number of samples, etc.)

        .. code-block::

            ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            sample                         |    mutation_mode   | min_sample_reads   | # of cohort samples   |   reads per cohort sample  |   # of normal samples allowed   |   normal_cohort_id   |    reads per normal sample   |   init_cancer   |   filter_sample   |   filter_sample_cohort   |   filter_sample_cohort_cohortnormal   |   filter_sample_cohort_cohortnormal_uniprot   |   info              |
            ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            TCGA-OR-A5J1-01A-11R-A29R-07   |    somatic         |      3.0           |     10                |      0.0                   |             1                   |         GtexCore     |            0.0               |      1865044    |      1661806      |      1661432             |                  21399                |                    21378                      |    GtexcoreCore     |
            TCGA-A2-BCRT-01A-11R-A084-07   |       ref          |      2.0           |      1                |      0.0                   |             2                   |         GtexCore     |            10.0              |      1865044    |      1661806      |      1661432             |                  21401                |                    21380                      |    GtexcoreCore     |
            ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


.. _mhcbind_out:

Outputs mode `mhcbind`
----------------------

1. **[--input-peptides-file]**: Intermediate file created if *--partitioned-tsv* is provided. It contains the unique set of kmers for which the user wants to predict binding to the MHC complex.

    **Format:** It will be a csv file, and it will contain a column with the unique kmers, without a header or index.

    **Pipeline relevance:** It contains the unique set of kmers present in the partitioned tsv files provided under --partitioned-tsv. This corresponds to the files 1 and 2 found in the :ref:`output section <output-tsv-cancerspecif>`.

    **Technical note:** The intermediate file will be stored in the path provided under *--input-peptides-file* in *--argstring*. This file will be passed to the mhcbind tool argument *--input-peptides-file*.

    .. code-block::

        AAGDDENHN
        AAPGQHLQA
        ACSNFIFKH


2. **[--output-csv]**:  Output file generated by the selected MHC tool.

    **Format:** csv file, and the exact format will depend on the MHC tool selected.

    **Technical note:** The path provided as *--output-csv* by the user in the *--argstring* will be passed to the mhcbind tool argument *--output-csv*.

    **Note**: The following example has the format of mhcflurry tool.

    .. code-block::

       +----------------------+--------+----------+-------------+---------------------+--------------------+-------------------+------------------------+-------+
       | source_sequence_name | offset | peptide  | allele      | score               | affinity           |percentile_rank    | prediction_method_name |length |
       +----------------------+--------+----------+-------------+---------------------+--------------------+-------------------+------------------------+-------+
       | Na                   | 0      | AAGDDENH | HLA-A*02:01 | 0.078546            | 30000.4500000      | 37.000003         | mhcflurry              | 9     |
       | Na                   | 0      | AAPGQHLQ | HLA-A*02:01 | 0.314567            | 900.89675          | 2.08654           | mhcflurry              | 9     |
       | Na                   | 0      | ACSNFIFK | HLA-A*02:01 | 0.000044            | 234.567345         | 50.05467392       | mhcflurry              | 9     |
       +----------------------+--------+----------+-------------+---------------------+--------------------+-------------------+------------------------+-------+

3. **[--output-csv]_With{--bind-score-method}LessLim{--bind-score-threshold}.tsv**: File containing the filtered version of **[--output-csv]**. It is generated if filtering options are provided, and *--less-than* is set to True. It will contain kmers that have a score *--bind-score-method* smaller than *--bind-score-threshold*.

    .. code-block::

       +----------------------+--------+----------+-------------+---------------------+--------------------+-------------------+------------------------+-------+
       | source_sequence_name | offset | peptide  | allele      | score               | affinity           |percentile_rank    | prediction_method_name |length |
       +----------------------+--------+----------+-------------+---------------------+--------------------+-------------------+------------------------+-------+
       | Na                   | 0      | AAGDDENH | HLA-A*02:01 | 0.078546            | 30000.4500000      | 37.000003         | mhcflurry              | 9     |
       +----------------------+--------+----------+-------------+---------------------+--------------------+-------------------+------------------------+-------+

4. **[--output-csv]_With{--bind-score-method}MoreLim{--bind-score-threshold}.tsv**: File containing the filtered version of **[--output-csv]**. It is generated if filtering options are provided, and *--less-than* is set to False. It will contain kmers that have a score *--bind-score-method* bigger than *--bind-score-threshold*. **Format:** Same format as the file **[--output-csv]_With{--bind-score-method}LessLim{--bind-score-threshold}.tsv**

.. _pepquery_out:

Outputs mode `pepquery`
------------------------

1. **[-i]**: This file is created if *--partitioned-tsv* is provided. It contains the expanded peptides generated from the kmers obtained in the cancerspecif mode. This will be the file introduced as an input to the PepQuery validation tool. The name and saving path will correspond with the information introduced under -i in the --argstring. **Format**: File with a peptide per row, no header.

    .. code-block::

            LSLVHPGTRRITKRRRQYPYVIASCQREAGCRGIICS
            NEVIGECIACSASFDATTIGRSRHRESSSLVSDGWACRGSATARPPNPRRAVLCKSIEPTYG
            YPYVIASCQREAGCRGIICS

2. **[-o]**: This is a folder containing the raw results of the pepQuery software tool. If the user is interested in accessing the results directly generated from the software, one can look in this folder. The name and saving directory correspond to the directory introduced under -o in the --argstring. A more detailed description of the contents can be found in the `PepQuery documentation <http://pepquery.org/document.html#saoutput>`_.

3. **[--output_dir]/peptides_validated.tsv.gz**: This file contains the main results of pepQuery formatted in a more interpretable way. In the file, the user will encounter different columns:

    - **peptide**: The peptide sequence under study. If --partitioned-tsv was provided this peptide sequence corresponds to the expanded kmer sequence from cancerspecif output kmers.
    - **modification**: The modification present on the peptide.
    - **spectrum**: MS/MS spectrum that the peptide matched with a sufficiently high score.
    - **score**: The score of the peptide-spectrum match.
    - **confident**: Whether the match is confident or not. In order for a match to be considered confident it has to pass the filtering steps described in the modes section.
    - **pvalue**: pvalue of the statistical evaluation step.
    - **# of reference DB peptides matching spectrum better than study peptide**: Number of peptides in the reference database that match the spectrum better than the study peptide. If there is one reference peptide matching better, the match will be considered non-confident. The following criteria will not be assessed and a NaN will appear in the next columns.
    - **# random shuffled peptides matching spectrum better that study peptide**:Number of random shuffled peptides matching the spectrum better than the study peptide. If there is a number of random shuffled peptides that match the spectra better, such that the pvalue > 0.01, the match will be considered non-confident. The following criteria will not be assessed and a NaN will appear in the next columns.
    - **# ptm-modified proteins matching better the spectra filtering summary that study peptide**: Number of post translational modified proteins of the reference databse that match better the spectrum than the study peptide. If there is one post translationally modified reference peptide matching better, the match will be considered non-confident.
    - **filtering summary**: Summary of the filtering steps. It will show whether the peptide passed all the filters or in which specific filtering step it failed.

    .. code-block::

       +--------------------------------------------+------------------------------------------------------------------------------------------------------------------+-------------------------------------------+--------------------+---------------------+--------------+-------------------------------------------------------------------------+------------------------------------------------------------------------+--------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+
       | peptide                                    | modification                                                                                                     | spectrum                                  | score              | confident           | pvalue       |# of reference DB peptides matching spectrum better than study peptide   | # random shuffled peptides matching spectrum better that study peptide |  # ptm-modified proteins matching better the spectra filtering summary that study peptide  |      filtering summary                                                                      |
       +--------------------------------------------+------------------------------------------------------------------------------------------------------------------+-------------------------------------------+--------------------+---------------------+--------------+-------------------------------------------------------------------------+------------------------------------------------------------------------+--------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+
       | VPEPGCTKVPEPGCTKFPEPGYTKVPVPGYTKVPEPCPSTVT | Carbamidomethylation of C@6[57.0215];Carbamidomethylation of C@14[57.0215];Carbamidomethylation of C@37[57.0215] | 03CPTAC_LUAD_W_BI_20180521_KR_f13:25749:2 |  25.64554920977335 |  Yes                | 0.005        |                       0                                                 |                              20                                        |                        0                                                                   | The peptide passed all the filters and the identified spectra is considered confident       |
       | LELMSVLSSGSLVHSRSSMLRRDH                   | Oxidation of M@4[15.9949];Oxidation of M@19[15.9949]                                                             | 03CPTAC_LUAD_W_BI_20180521_KR_f13:25749:2 |  14.44877947689417 |  No                 | 0.5566       |                       1                                                 |                              NaN                                       |                        NaN                                                                 | Failed at competitive filtering based on reference sequences (step 3).                      |
       | MSSQQQKQPCIPPPQLQQQQVKQPCQPPPQ             | Carbamidomethylation of C@13[57.0215]                                                                            | 03CPTAC_LUAD_W_BI_20180521_KR_f13:25749:2 |  12.510356427561199|  No                 | 0.12333      |                       0                                                 |                              146                                       |                        NaN                                                                 | Failed at the statistical evaluation based on random shuffling (step 4). The pvalue is >0.01|
       +--------------------------------------------+------------------------------------------------------------------------------------------------------------------+-------------------------------------------+--------------------+---------------------+--------------+-------------------------------------------------------------------------+------------------------------------------------------------------------+--------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+


