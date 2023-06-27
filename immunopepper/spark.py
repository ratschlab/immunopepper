import logging
import numpy as np
import os
import pandas as pd
import pathlib
import sys
from pyspark.sql import functions as sf
from pyspark.sql import types as st


def process_libsize(path_lib, custom_normalizer):
    '''
    Loads and returns the normalisation values per sample
    Normalisation formulae: (count / 75 quantile expression in sample) * A
    If no normalisation factor is provided: A = median across samples
    If normalisation factor is provided: A = normalisation factor
    :param path_lib: str path with library size file
    :param custom_normalizer: custum normalisation factor
    :return: dataframe with 75 quantile expression in sample / A values
    '''
    lib = pd.read_csv(path_lib, sep='\t')
    if custom_normalizer:
        lib['libsize_75percent'] = lib['libsize_75percent'] / custom_normalizer
    else:
        lib['libsize_75percent'] = lib['libsize_75percent'] / np.median(lib['libsize_75percent'])
    lib['sample'] = [sample.replace('-', '').replace('.', '').replace('_','') for sample in lib['sample']]
    lib = lib.set_index('sample')
    return lib


def process_build_outputs(spark, index_name, jct_col, jct_annot_col, rf_annot_col,
                        path_matrix=None, whitelist=None, cross_junction=None,
                        filterNeojuncCoord=None, filterAnnotatedRF=None,
                        output_dir=None, separate_back_annot=None, tot_batches=None, batch_id=None):
    '''
    Various preprocessing steps including
    - Loading of data from immunopepper mode build
    - Separation of the kmers derived solely from the annotation (expression is null in samples)
    - Filtering on junction and reading frame annotated status
    - Select whitelist samples
    :param spark: spark context
    :param index_name: str kmer column name
    :param jct_col: str junction column name
    :param jct_annot_col: str junction annotated flag column name
    :param rf_annot_col: str reading frame annotated flag column name
    :param path_matrix: str path for multi-sample count matrix from immunopepper mode build:  kmer|is_cross_junction|junctionAnnotated|readFrameAnnotated|sample_1|sample_2|...
    :param path_kmer_file: str path for single sample kmer file from immunopepper mode build:  kmer|id|segmentexpr|iscrossjunction|junctionexpr|junctionAnnotated|readFrameAnnotated|
    :param col_expr_kmer_file: str name of the expression column of interest in path_kmer_file. e.g. 'segmentexpr' or 'junctionexpr'
    :param target_sample: str sample for which the path_kmer_file has been computed
    :param whitelist: list whitelist for samples
    :param cross_junction: bool whether the count matrix contains junction counts or segment counts
    :param annot_flag: list with the intruction codes on how to treat the reading frame and junction annotated flags
    :param output_dir: str path for output directory
    :param filterNeojuncCoord: bool if True, filter for kmers from neojunctions,
     i.e. with non-annotated junction coordinates
    :param filterAnnotatedRF: bool if True, filter for kmers from annotated reading frames,
     i.e. with reading frames found in transcripts from the annotation, not propagated
    :param separate_back_annot: str, if provided, the kmers derived from the annotation only (i.e. without read support
    in any cohort sample) are saved into this path.
    - The kmers are thereby skipped from expression threshold filtering in the background matrix.
    - They will be by default removed from the foreground matrix as "annotated kmers"
    :param tot_batches: int batch mode only, total number of batches
    :param batch_id: int batch mode only, id of the batch
    :return: Preprocessed spark dataframe
    '''

    def RenamedCols(df):
        old_name = df.columns
        new_names = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in old_name]
        df = df.toDF(*new_names)
        return df


    def filter_whitelist(matrix, name_list, index_name, jct_annot_col, rf_annot_col):
        name_list.extend([index_name, jct_annot_col, rf_annot_col])
        return matrix.select([sf.col(name_) for name_ in name_list])


    if (path_matrix is not None):
        # Load immunopepper kmer candidates
        rename = False  # For development
        logging.info(f'Load input {path_matrix}')
        matrix = loader(spark, path_matrix, header=True)
        partitions_ = matrix.rdd.getNumPartitions()
        logging.info(f'...partitions: {partitions_}')

        if rename:
            logging.info("Rename")
            matrix = RenamedCols(matrix)

        # In batch mode: Dev remove kmers for the matrix based on the hash function
        if tot_batches:
            logging.info((f'Filter foreground and background based on Hash ; '
                          f'making {tot_batches} batches, select batch number {batch_id}'))

            matrix = matrix.filter(sf.abs(sf.hash('kmer') % tot_batches) == batch_id)

        # Filter on junction status depending on the content on the matrix
        matrix = matrix.drop(jct_col)

        # Separate kmers only present in the backbone annotation from the ones supported by the reads of any sample
        partitions_ = matrix.rdd.getNumPartitions()
        logging.info(f'...partitions: {partitions_}')
        if separate_back_annot:
            logging.info("Isolating kmers only in backbone annotation")
            matrix = split_only_found_annotation_backbone(separate_back_annot, output_dir, matrix, index_name,
                                                          jct_annot_col, rf_annot_col)

        # Filter according to annotation flag
        matrix = filter_on_junction_kmer_annotated_flag(matrix, jct_annot_col, rf_annot_col,
                                                      filterNeojuncCoord, filterAnnotatedRF)

        # Reduce samples (columns) to whitelist
        logging.info("Cast types")
        if whitelist is not None:
            whitelist = pd.read_csv(whitelist, sep='\t', header=None)[0].to_list()
            whitelist = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in whitelist]
            matrix = filter_whitelist(matrix, whitelist, index_name, jct_annot_col, rf_annot_col)

        partitions_ = matrix.rdd.getNumPartitions()
        logging.info(f'...partitions: {partitions_}')
        return matrix
    else:
        return None


def split_only_found_annotation_backbone(separate_back_annot, output_dir, matrix, index_name,
                                         jct_annot_col, rf_annot_col):
    '''
    Separate kmers only present in the backbone annotation from the ones supported by the reads of any sample:
    A kmer is solely derived from the backbone annotation if (all samples have zero expression)
    The opposite condition is implemented to make use of short-circuit evaluation

    :param matrix: spark dataframe matrix to filter
    :param index_name: str kmer column name
    :param jct_annot_col: string junction is annotated column
    :param rf_annot_col: string reading frame is annotated column
    :return: spark matrix with expressed kmers, spark serie with kmers only in backbone annotation
    '''
    expressed = ' OR '.join( [f'({col_name} != 0.0)'
                                   for col_name in matrix.schema.names
                                   if col_name not in [index_name, jct_annot_col, rf_annot_col]]) # SQL style because many cols
    matrix = matrix.withColumn('expressed',
                             sf.when(sf.expr(f"{expressed}"), True).otherwise(False))

    kmers_AnnotOnly = matrix.select(index_name).where(sf.col('expressed') == False)
    matrix_Expressed = matrix.filter(sf.col('expressed') == True)
    matrix_Expressed = matrix_Expressed.drop(sf.col('expressed'))
    save_spark(kmers_AnnotOnly, output_dir, separate_back_annot)
    return matrix_Expressed


def filter_on_junction_kmer_annotated_flag(matrix, jct_annot_col, rf_annot_col, filterNeojuncCoord, filterAnnotatedRF):
    '''
    Filters according to the junction and kmer annotated flag
    :param matrix: spark dataframe matrix to filter
    :param jct_annot_col: string junction is annotated column
    :param rf_annot_col: string reading frame is annotated column
    :param filterNeojuncCoord: bool if True, filter for kmers from neojunctions,
     i.e. with non-annotated junction coordinates
    :param filterAnnotatedRF: bool if True, filter for kmers from annotated reading frames,
     i.e. with reading frames found in transcripts from the annotation, not propagated
    '''
    # Keep k-mers according to annotation flag
    if filterNeojuncCoord:
        matrix = matrix.filter(f'({jct_annot_col} == {False})')
    if filterAnnotatedRF:
        matrix = matrix.filter(f'({rf_annot_col} == {True})')
    return matrix


def combine_normals(normal_segm, normal_junc):
    ''' Given two matrices [kmers] x [samples, metadata] containing segment expression and junction expression,
    takes the union expression of the two matrices
    :param normal_segm: spark dataframe matrix with segment expression counts
    :param normal_junc: spark dataframe matrix with junction expression counts
    :return: spark dataframe expression matrix after combining junction and segment expression
    '''
    if (normal_segm is not None) and (normal_junc is not None):
        logging.info("Combine segment and edges")
        normal_matrix = normal_segm.union(normal_junc)
        return normal_matrix
    elif normal_segm is None:
        return normal_junc
    elif normal_junc is None:
        return normal_segm


def check_interm_files(out_dir, expr_limit, n_samples_lim, target_sample='', tag='normals', batch_tag=''):
    '''
    Filtering steps for normal (resp. cancer) samples are saved as intermediate files because it is an expensive operation
    The function checks the presence of the intermediate filtering files to decide whether to perform
    - the full preprocessing + threshold filtering steps
    - or simply re-load the intermediate files
    :param out_dir: str path for output directory
    :param expr_limit: float expression limit threshold to keep a kmer
    :param n_samples_lim: int number of samples that need to pass the expression limit
    :param target_sample: str name of the sample of interest.
    To be excluded in the number of samples that pass the expression limit
    :param tag: str tag related to the type of samples. Example cancer or normal
    :param batch_tag: str batch mode, batch tag to be appended to intermediate file
    :returns:
    - launch_preprocess: bool, whether to perform the full preprocessing + threshold filtering steps
     or simply re-load the intermediate files
    - path_interm_matrix_for_express_threshold, path_interm_matrix_for_sample_threshold,
     path_interm_kmers_annotOnly are respectively the path (str) where
     the expression-filtered matrix, the sample-filtered matrix and
     the kmers derived solely from the annotation are saved
    '''
    base_n_samples = 1
    base_expr = 0.0
    format_tag = '.tsv.gz'
    # For cancer matrix intermediate file the recurrence filter is not applied to the target sample
    if target_sample:
        suffix = f'Except{target_sample}'
    else:
        suffix = ''

    # For normal samples the expression threshold filtering is not applied to the kmers found only in the annotation
    # but not in the background samples. These kmers will be by default removed from the foreground matrix.
    if tag == 'normals':
        path_interm_kmers_annotOnly = os.path.join(out_dir, f'kmers_derived_solely_from_annotation{format_tag}')
    else:
        path_interm_kmers_annotOnly = None

    # Saving paths

    path_interm_matrix_for_express_threshold = os.path.join(out_dir,
                          f'interm_{tag}_combiExprCohortLim{expr_limit}Across{base_n_samples}{suffix}{batch_tag}{format_tag}')
    path_interm_matrix_for_sample_threshold = os.path.join(out_dir,
                              f'interm_{tag}_combiExprCohortLim{base_expr}Across{base_n_samples}{suffix}{batch_tag}{format_tag}')
    # Check existence
    if (expr_limit and os.path.isfile(os.path.join(path_interm_matrix_for_express_threshold, '_SUCCESS'))) \
        and (n_samples_lim is not None and os.path.isfile(os.path.join(path_interm_matrix_for_sample_threshold, '_SUCCESS'))):

        logging.info((f'Intermediate {tag} filtering already performed in: {path_interm_matrix_for_express_threshold} '
                      f' and {path_interm_matrix_for_sample_threshold}. Re-loading {tag} intermediate data...'))
        logging.info((f'Proceed with care! Using intermediate files means ignoring --filterNeojuncCoord, '
                      f'--filterAnnotatedRF parameter.'))
        launch_preprocess = False
    else:
        logging.info(f'At least one intermediate {tag} filtering file is missing.')
        logging.info(f'Will compute full filtering steps according to user input parameters')
        launch_preprocess = True

    return launch_preprocess, \
           path_interm_matrix_for_express_threshold, \
           path_interm_matrix_for_sample_threshold, \
           path_interm_kmers_annotOnly


def filter_hard_threshold(matrix, index_name, jct_annot_col, rf_annot_col, libsize,
                          expr_limit, n_samples_lim, path_e, path_s, target_sample='',
                          tag='normals', on_the_fly=False):
    '''
    Filter samples based on >0 and >=X reads expressed (expensive operations) and save intermediate files.
    Additional thresholds will be applied specifically for cancer or normals matrices in subsequent combine functions
    a. kmer need to be if >= X reads in >= 1 sample -> saved as intermediate file as path_e (expression)
    b. kmer needs to be >0 reads in >= 1 sample -> saved as intermediate file as path_s (sample)
    The expressions are normalized for library size.
    The filtering is either performed on the full cohort in the matrix or in the cohort excluding a target sample.
    :param matrix: spark dataframe for count matrix
    :param index_name: str kmer column name
    :param jct_annot_col: str junction annotated flag column name
    :param rf_annot_col: str reading frame annotated flag column name
    :param libsize: dataframe with library size
    :param expr_limit: float expression limit threshold to keep a kmer
    :param n_samples_lim: int number of samples that need to pass the expression limit
    :param path_e: path where the expression-filtered matrix (a) is saved
    :param path_s: path where the sample-filtered matrix (b) is saved
    :param target_sample: str name of the sample of interest.
     To be excluded in the number of samples that pass the expression limit
    :param tag: str tag related to the type of samples. Example cancer or normal
    :param on_the_fly: bool. whether to save intermediate file after counting the number of columns passing a threshold
    :returns
        matrix_e: RDD. intermediate object after applying >= X reads in >= 1 sample
        matrix_s: RDD. intermediate object after applying >0 reads in >= 1 sample
    '''

    base_n_samples = 1
    base_expr = 0.0
    matrix_e, matrix_s = None, None


    if target_sample:
        logging.info(f'Target sample {target_sample} not included in the cohort filtering')
    if libsize is not None:
        matrix = matrix.select(index_name, *[
            sf.round(sf.col(name_) / libsize.loc[name_, "libsize_75percent"], 2).alias(name_)
            for name_ in matrix.schema.names if name_ not in [index_name, 'coord', jct_annot_col, rf_annot_col]])

    # Expression filtering, take k-mers with >= X reads in >= 1 sample
    if expr_limit and (not os.path.isfile(os.path.join(path_e, '_SUCCESS'))):
        logging.info(f'Filter matrix with cohort expression support >= {expr_limit} in {base_n_samples} sample')
        # Fill the expression matrix with 1 if expression threshold is met, 0 otherwise
        # Skip target sample and metadata
        matrix_e = matrix.select(index_name, *[
            sf.when(sf.col(name_) >= expr_limit, 1).otherwise(0).alias(name_)
            for name_ in matrix.schema.names if name_
            not in [target_sample, index_name, jct_annot_col, rf_annot_col]])
        # Map each kmer: x[0] to the number of samples where the expression threshold is met: sum(x[1:])
        # Get a tuple (kmer, number of samples where expression >= threshold)
        # Then filter kmers based on the number of samples
        # This uses rdds syntax because spark dataframe operations are slower
        matrix_e = matrix_e.rdd.map(tuple).map(lambda x: (x[0], sum(x[1:]))).\
            filter(lambda x: x[1] >= base_n_samples)
        if not on_the_fly:
            logging.info(f'Save intermediate 1/2 {tag} filtering file to {path_e}')
            matrix_e.map(lambda x: "%s\t%s" % (x[0], x[1])).saveAsTextFile(path_e, \
                          compressionCodecClass="org.apache.hadoop.io.compress.GzipCodec")

    # Sample filtering, take k-mers with exclude >0 reads in >= 1 sample
    if (n_samples_lim is not None) and (not os.path.isfile(os.path.join(path_s, '_SUCCESS'))):
        logging.info(f'Filter matrix with cohort expression support > {base_expr} in {base_n_samples} sample')
        # Fill the expression matrix with 1 if expression threshold is met, 0 otherwise
        # Skip target sample and metadata
        matrix_s = matrix.select(index_name, *[
            sf.when(sf.col(name_) > base_expr, 1).otherwise(0).alias(name_)
            for name_ in matrix.schema.names if name_
            not in [target_sample, index_name, jct_annot_col, rf_annot_col]])
        # Map each kmer: x[0] to the number of samples where the expression threshold is met: sum(x[1:])
        # Get a tuple (kmer, number of samples where expression >= threshold)
        # Then filter kmers based on the number of samples
        # This uses rdds syntax because spark dataframe operations are slower
        matrix_s = matrix_s.rdd.map(tuple).map(lambda x: (x[0], sum(x[1:]))).\
            filter(lambda x: x[1] >= base_n_samples)
        if not on_the_fly:
            logging.info(f'Save intermediate 2/2 {tag} filtering file to {path_s}')
            matrix_s.map(lambda x: "%s\t%s" % (x[0], x[1])).saveAsTextFile(path_s,  \
                          compressionCodecClass="org.apache.hadoop.io.compress.GzipCodec")

    return matrix_e, matrix_s




def combine_hard_threshold_normals(spark, path_normal_kmers_e, path_normal_kmers_s,
                                   normal_matrix_e, normal_matrix_s,
                                   n_samples_lim_normal, index_name):
    '''
    Filter samples based on X reads and H samples.
    a. kmer need to be if >= X reads in >= 1 sample -> load path_*_e (expression) file
    b. kmer needs to be >0 reads in >= H samples -> load path_*_s (sample) and apply H threshold on the fly
    The kmers selected need to pass a. OR b.
    :param spark: spark context
    :param path_normal_kmers_e: str path for intermediate file which applied >= X reads in >= 1 sample
    :param path_normal_kmers_s: str path for intermediate file which applied >0 reads in >= 1 sample
    :param normal_matrix_e: RDD. intermediate object after applying >= X reads in >= 1 sample
    :param normal_matrix_s: RDD. intermediate object after applying >0 reads in >= 1 sample
    :param n_samples_lim_normal: int number of samples in which any number of reads is required
    :param index_name: index_name: str kmer column name
    :return: spark dataframe filtered with >= X reads in 1 sample OR >0 reads in H samples
    '''
    #  Convert or re-load matrix expression threshold (Counts >= X reads in >= 1 sample)
    if normal_matrix_e:
        normal_matrix_e = normal_matrix_e.toDF(['kmer', 'n_samples'])
        normal_matrix_e = normal_matrix_e.select(sf.col(index_name))
    elif path_normal_kmers_e:
        if (path_normal_kmers_s != path_normal_kmers_e):
            normal_matrix_e = spark.read.csv(path_normal_kmers_e, sep=r'\t', header=False)
            normal_matrix_e = normal_matrix_e.withColumnRenamed('_c0', index_name)
            normal_matrix_e = normal_matrix_e.select(sf.col(index_name))


    # convert or re-load matrix sample threshold (Counts >0  reads in >= 1 sample)
    if normal_matrix_s:
        normal_matrix_s = normal_matrix_s.toDF(['kmer', 'n_samples'])
    elif path_normal_kmers_s:
        normal_matrix_s = spark.read.csv(path_normal_kmers_s, sep=r'\t', header=False)
        normal_matrix_s = normal_matrix_s.withColumnRenamed('_c0', index_name)
        normal_matrix_s = normal_matrix_s.withColumnRenamed('_c1', "n_samples")


    # Threshold on number of samples
    if normal_matrix_s:
        if n_samples_lim_normal > 1: #(Counts >0  reads in >= H samples)
            logging.info(f'Filter matrix with cohort expression support > {0} in {n_samples_lim_normal} sample(s)')
            normal_matrix_s = normal_matrix_s.filter(sf.col('n_samples') >= n_samples_lim_normal)
        normal_matrix_s = normal_matrix_s.select(sf.col(index_name))

    if normal_matrix_e and not normal_matrix_s:
        return normal_matrix_e
    elif normal_matrix_s and not normal_matrix_e:
        return normal_matrix_s
    elif normal_matrix_e and normal_matrix_s:
        normal_matrix_res = normal_matrix_e.union(normal_matrix_s) # Do not make distinct
        return normal_matrix_res
    else:
        return None


def combine_hard_threshold_cancers(spark, cancer_matrix, path_cancer_kmers_e, path_cancer_kmers_s,
                                   inter_matrix_expr_c, inter_matrix_sample_c,
                                   cohort_expr_support_cancer, n_samples_lim_cancer, index_name):
    '''
     Filter samples based on X reads and H samples.
    a. kmer need to be expressed with >= X reads and this in >= H samples
    Will be performed by loading path_*_e (expression) file or path_*_s (sample) and applying H threshold on the fly
    Then set a. will be applied to preprocessed matrix
    :param spark: spark context
    :param cancer_matrix: spark dataframe with preprocessed foreground
    :param path_normal_kmers_e: str path for intermediate file which applied >= X reads in >= 1 sample
    :param path_normal_kmers_s: str path for intermediate file which applied >0 reads in >= 1 sample
    :param inter_matrix_expr_c: RDD. intermediate object after applying >= X reads in >= 1 sample
    :param inter_matrix_sample_c: RDD. intermediate object after applying >0 reads in >= 1 sample
    :param cohort_expr_support_cancer: float expression threshold to be met for the kmer
    :param n_samples_lim_cancer: int number of samples in which the expression threshold need to be met
    :param index_name: str kmer column name
    :return: spark dataframe with foreground filtered for >= X reads and in >= H samples
    '''
    # Convert or re-load intermediate files
    if cohort_expr_support_cancer != 0.0:
        if inter_matrix_expr_c:
            valid_foreground = inter_matrix_expr_c.toDF(['kmer', 'n_samples'])
        elif path_cancer_kmers_e:
            valid_foreground = spark.read.csv(path_cancer_kmers_e, sep=r'\t', header=False)
            valid_foreground = valid_foreground.withColumnRenamed('_c0', index_name).withColumnRenamed('_c1', "n_samples")

    else:
        if inter_matrix_sample_c:
            valid_foreground = inter_matrix_sample_c.toDF(['kmer', 'n_samples'])
        elif path_cancer_kmers_s:
            valid_foreground = spark.read.csv(path_cancer_kmers_s, sep=r'\t', header=False)
            valid_foreground = valid_foreground.withColumnRenamed('_c0', index_name).withColumnRenamed('_c1', "n_samples")


    # kmer need to be expressed with >= X reads and this in >= H samples
    if n_samples_lim_cancer > 1:
        logging.info( (f'Filter matrix with cohort expression support >= {cohort_expr_support_cancer} '
                       f'in {n_samples_lim_cancer} sample(s)'))
        valid_foreground = valid_foreground.filter(sf.col('n_samples') >= n_samples_lim_cancer)
    valid_foreground = valid_foreground.select(sf.col(index_name))

    # Apply to preprocessed matrix
    cancer_cross_filter = cancer_matrix.join(valid_foreground, ["kmer"],
                                             how='right')
    return cancer_cross_filter


def filter_expr_kmer(matrix_kmers, filter_field, threshold, libsize=None):
    '''
    Filters a spark dataframe on a threshold with or without prior normalization
    :param matrix_kmers: spark dataframe with expression counts
    :param filter_field: str name of column on which to filter
    :param threshold: float threshold (>= threshold)
    :param libsize_c: dataframe with library size
    :return: filtered spark dataframe with expression counts
    '''
    partitions_ = matrix_kmers.rdd.getNumPartitions()
    logging.info(f'...partitions: {partitions_}')
    if threshold != 0:
        logging.info(f'Filter with {filter_field} >= {threshold}')
        # Normalize by library size
        if libsize is not None:
            logging.info("Normalizing cancer counts")
            matrix_kmers = matrix_kmers.filter(sf.round(sf.col(filter_field)/
                                                        libsize.loc[filter_field, "libsize_75percent"], 2) >= threshold)
        else:
            matrix_kmers = matrix_kmers.filter(sf.col(filter_field) >= threshold)

    else:
        logging.info(f'Filter with {filter_field} > {threshold}')
        matrix_kmers = matrix_kmers.filter(sf.col(filter_field) > threshold)
    return matrix_kmers


def combine_cancer(cancer_kmers_segm, cancer_kmers_edge, index_name):
    '''
    Given two matrices [kmers] x [samples, metadata] containing segment expression and junction expression,
    takes the junction expression if available, otherwise take the segment expression
    (Note: Junction expression is the best expression proxy for junction kmers,
    therefore it is advised to provide only a junction expression matrix so that the software skips this step)
    :param cancer_kmers_segm: spark dataframe matrix with segment expression counts
    :param cancer_kmers_edge: spark dataframe matrix with junction expression counts
    :param index_name: str name of the kmer column
    :return: spark dataframe matrix with combined expression counts
    '''
    if (cancer_kmers_segm is not None) and (cancer_kmers_edge is not None):
        logging.info("Combine segment and edges")
        cancer_kmers_segm = cancer_kmers_segm.join(cancer_kmers_edge,
                                                   cancer_kmers_segm[index_name] == cancer_kmers_edge[index_name],
                                                   how='left_anti')
        cancer_kmers = cancer_kmers_edge.union(cancer_kmers_segm)
        partitions_ = cancer_kmers.rdd.getNumPartitions()
        logging.info(f'...partitions cancer filtered: {partitions_}')
        return cancer_kmers
    elif (cancer_kmers_segm is None):
        return cancer_kmers_edge
    elif (cancer_kmers_edge is None):
        return cancer_kmers_segm


def remove_external_kmer_list(spark, path_external_kmer_database, spark_matrix, index_name, header=False):
    '''
    Takes an expression matrix of the format [kmers] x [samples, metadata] and removes an external database
    :param spark: spark context
    :param path_normal_kmer_list: str. path of file containing kmers to be removed
    :param spark_matrix: str path for count matrix
    :param index_name: str kmer column name
    :param header: bool. whether the file contains a header
    :return: Filtered matrix for external kmer database
    '''
    logging.info("Load input {}".format(path_external_kmer_database))
    external_database = loader(spark, path_external_kmer_database, header=header)
    external_database = external_database.select(sf.col(index_name))
    if spark_matrix:
        spark_matrix = spark_matrix.union(external_database)
    else:
        spark_matrix = external_database
    return spark_matrix


def remove_uniprot(spark, cancer_kmers, uniprot, index_name):
    '''
    Filters a spark dataframe against uniprot or any peptide database.
    Equivalence between leucine and isoleucine assumed.
    :param spark: spark context
    :param cancer_kmers: spark dataframe matrix with expression counts for cancer
    :param uniprot: str path for uniprot file
    :param index_name: str name of the kmer column
    :return: spark dataframe matrix for cancer after uniprot filtering
    '''
    def I_L_replace(value):
        return value.replace('I', 'L')
    if uniprot is not None:
        logging.info("Filter out uniprot")
        uniprot = spark.read.csv(uniprot, sep='\t', header=None)
        uniprot_header = index_name + "_IL_eq"
        uniprot = uniprot.withColumnRenamed("_c0", uniprot_header)
        # Make isoleucine and leucine equivalent
        I_L_equ = sf.udf(I_L_replace, st.StringType())
        uniprot = uniprot.withColumn(uniprot_header, I_L_equ(uniprot_header))
        cancer_kmers = cancer_kmers.withColumn(uniprot_header, I_L_equ(index_name))
        cancer_kmers = cancer_kmers.join(uniprot, cancer_kmers[uniprot_header] == uniprot[uniprot_header],
                                         how='left_anti')
    return cancer_kmers


def save_spark(cancer_kmers, output_dir, path_final_fil, outpartitions=None):
    '''
    Saves a spark dataframe as a single or partitioned csv file
    :param cancer_kmers: spark dataframe matrix with expression counts for cancer
    :param output_dir: str path for output directory
    :param path_final_fil: str path to save the spark dataframe
    :param outpartitions: int number of partitions for saving
    '''
    # save
    logging.info(f'>>>> Save to {path_final_fil}')
    pathlib.Path(output_dir).mkdir(exist_ok=True, parents=True)
    if outpartitions is not None:
        cancer_kmers.repartition(outpartitions).write.mode('overwrite')\
            .options(header="true", sep="\t", compression="gzip").format("tsv.gz").csv(path_final_fil)
    else:
        cancer_kmers.write.mode('overwrite')\
            .options(header="true", sep="\t", compression="gzip").format("tsv.gz").csv(path_final_fil)


def loader(spark, path_kmer, header=False):
    '''
    Loads a parquet, csv or tsv (partitioned) file as a spark dataframe
    The user can provide either
    - a list of files with the same schema
    - a string. folder containing the files
    :param spark: spark context
    :param path_kmer: Either (1) a list of files or (b) a string with the basefolder where the files are placed
    :param header: bool. whether the file contains a header
    :return: loaded spark dataframe

    '''
    if ('parquet' in path_kmer) or ('pq' in path_kmer): #single file or folder
        matrix = spark.read.parquet(path_kmer)
    elif ('parquet' in path_kmer[0]) or ('pq' in path_kmer[0]): # list
        matrix = spark.read.parquet(*path_kmer)
    elif ('csv' in path_kmer[0]) or ('csv' in path_kmer):
        matrix = spark.read.csv(path_kmer, sep=',', header=header)
    elif ('tsv' in path_kmer[0]) or ('tsv' in path_kmer):
        matrix = spark.read.csv(path_kmer, sep=r'\t', header=header)
    elif ('gzip' in path_kmer[0]) or ('gzip' in path_kmer): #If tsv or csv not in path -> assume gzip are tab separated
        matrix = spark.read.csv(path_kmer, sep=r'\t', header=header)
    elif ('gz' in path_kmer[0]) or ('gz' in path_kmer):
        matrix = spark.read.csv(path_kmer, sep=r'\t', header=header)
    else:
        logging.error(f'Cannot determine file type of {path_kmer}. Please include .parquet .pq .tsv .csv suffix to the files or the folder (partitionned)')
        sys.exit()


    return matrix


def output_count(perform_count, matrix, report_count, report_step, step_string):
    '''
    Performs a count operation on the number of kmers present in spark dataframe after a given filtering step
    Note: This operation is expensive but useful if the user is interested in intermediate filtering steps
    :param perform_count: bool whether to perform a count operation
    :param matrix: spark dataframe with kmer expression counts
    :param report_count: list to store result of successive counting operations
    :param report_step: list to store name of successive counting operations
    :param step_string: str name of the counting operation
    '''
    if perform_count:
        mycount = matrix.count()
        report_count.append(mycount)
        report_step.append(step_string)
        logging.info(f'# {step_string} n = {mycount} kmers')


def save_output_count(output_count, report_count, report_steps, prefix, cancer_sample_ori, mutation_mode,
                      sample_expr_support_cancer, cohort_expr_support_cancer, n_samples_lim_cancer,
                          cohort_expr_support_normal, n_samples_lim_normal, id_normals):
    '''
    Saves the number of kmers present in spark dataframe after each filtering step in a tabular file
    :param output_count: str path for count file of intermediate filtering steps
    :param report_count: list to store result of successive counting operations
    :param report_step: list to store name of successive counting operations
    :param prefix: str information to be added to the result line in an info column
    :param cancer_sample_ori: str id of target cancer sample which was filtered
    :param mutation_mode: str information about whether mutations where applied or not
    :param sample_expr_support_cancer: float normalized expression threshold for the cancer target sample
    :param cohort_expr_support_cancer: float normalized expression threshold for the cancer cohort
    excluding the target sample
    hich should be met in n samples
    :param n_samples_lim_cancer: int number of cancer samples in which the cancer cohort expression threshold
    should be met
    :param cohort_expr_support_normal: float normalized expression threshold for the normal cohort
    required in any sample (>=1)
    :param n_samples_lim_normal: int number of normal samples in which any number of reads is required (>0)
    :param id_normals: str id of the normal cohort (example gtex)
    '''
    if output_count:
        header = (f'{"sample"}\t{"mutation_mode"}\t{"min_sample_reads"}\t{"#_of_cohort_samples"}\t'
                  f'{"reads_per_cohort_sample"}\t{"#_normal_samples_allowed"}\t{"normal_cohort_id"}'
                  f'\t{"reads_per_normal_sample"}')
        line =   (f'{cancer_sample_ori}\t{mutation_mode}\t{sample_expr_support_cancer}\t{n_samples_lim_cancer}'
                  f'\t{cohort_expr_support_cancer}\t{n_samples_lim_normal}\t{id_normals}'
                  f'\t{cohort_expr_support_normal}')

        for idx in np.arange(len(report_count)):
            header += f'\t{report_steps[idx]}'
            line += f'\t{report_count[idx]}'
        if prefix:
            header += f'\t{"info"}'
            line += f'\t{prefix}'
        header += "\n"
        line += "\n"
        if not os.path.exists(output_count):
            with open(output_count,"w") as f:
                f.write(header)
        with open(output_count, "a") as f:
            f.write(line)
        logging.info(f'Save intermediate info to {output_count}')


def redirect_interm(interm_dir_norm, interm_dir_canc, output_dir):
    '''
    Set the directory to save intermediary file
    - The output directory
    - Any other specified normal or cancer directory
    Default. Uses output directory
    :param interm_dir_norm: str custom scatch dir path to save intermediate normal files
    :param interm_dir_canc: str custom scatch dir path to save intermediate cancer files
    :param output_dir: str output directory for the filtered matrix
    :return:
    '''

    if interm_dir_canc:
        cancer_out = interm_dir_canc
    else:
        cancer_out = output_dir
    if interm_dir_norm:
        normal_out = interm_dir_norm
    else:
        normal_out = output_dir
    return normal_out, cancer_out
