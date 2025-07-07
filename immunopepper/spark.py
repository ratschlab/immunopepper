import logging
import os
from immunopepper.sdisk import save_spark
from pyspark.sql import functions as sf
from pyspark.sql import types as st



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


def filter_hard_threshold(matrix, index_name, coord_name, jct_annot_col, rf_annot_col, libsize,
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
            for name_ in matrix.schema.names if name_ not in [index_name, coord_name, jct_annot_col, rf_annot_col]])

    # Expression filtering, take k-mers with >= X reads in >= 1 sample
    if expr_limit and (not os.path.isfile(os.path.join(path_e, '_SUCCESS'))):
        logging.info(f'Filter matrix with cohort expression support >= {expr_limit} in {base_n_samples} sample')
        # Fill the expression matrix with 1 if expression threshold is met, 0 otherwise
        # Skip target sample and metadata
        matrix_e = matrix.select(index_name, *[
            sf.when(sf.col(name_) >= expr_limit, 1).otherwise(0).alias(name_)
            for name_ in matrix.schema.names if name_
            not in [target_sample, index_name, coord_name, jct_annot_col, rf_annot_col]])
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
            not in [target_sample, index_name, coord_name, jct_annot_col, rf_annot_col]])
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





