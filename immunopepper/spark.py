import logging
import numpy as np
import os
import pandas as pd
import pathlib
from pyspark.sql import functions as sf
from pyspark.sql import types as st
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, Formula, r
from rpy2.robjects.packages import importr
import scipy
from scipy import stats

from immunopepper.io_ import save_pd_toparquet

pandas2ri.activate()
deseq = importr('DESeq2')
BiocParallel = importr('BiocParallel')
BiocGenerics = importr("BiocGenerics")



def DESeq2(count_matrix, design_matrix, normalize, cores=1): #TODO decide on statistical mode future
    '''
    Runs DEseq2
    :param count_matrix: dataframe matrix with counts
    :param design_matrix: dataframe  matrix with experiment design
    :param normalize: dataframe with library sizes
    :param cores: nt number of cores to use
    :return: dataframe with result of the negative binomial fit
    '''
    # gene_column = ''
    to_dataframe = ro.r('function(x) data.frame(x)')
    count_matrix = round(count_matrix)
    count_matrix = pandas2ri.py2rpy(count_matrix)
    design_matrix = pandas2ri.py2rpy(design_matrix)
    design_formula = Formula(' ~ 1')

    dds0 = deseq.DESeqDataSetFromMatrix(countData=count_matrix,
                                            colData=design_matrix,
                                            design=design_formula)
    dds0 = BiocGenerics.estimateSizeFactors(dds0, type="poscounts")
    order_size_factor = list(dds0.do_slot('colData').do_slot('rownames'))
    if normalize is not None:
        logging.info("Enforcing custom normalisation in DESeq2")
        dds0.do_slot('colData').do_slot('listData')[1] = ro.vectors.FloatVector(list(normalize.loc[order_size_factor,
                                                                                                   'libsize_75percent'])) # Enforce size factors
    else:
        logging.info("WARNING: default size factor of DESeq2 are used")
    dds = deseq.DESeq(dds0, parallel=True, BPPARAM=BiocParallel.MulticoreParam(cores),
                      sfType="poscounts", # Will run 1. estimation of size factors: estimateSizeFactors # parameter "poscounts"
                      fitType="parametric" # 2. estimation of dispersion: estimateDispersions # parameter "parametric"
                      )

    deseq_result = deseq.results(dds)
    fit_res = to_dataframe(deseq_result)
    disp = to_dataframe(deseq.dispersions(dds)).rename({'x': 'dispersion'}, axis = 1)
    disp.index = fit_res.index
    fit_res = pd.concat([fit_res['baseMean'], disp], axis=1)
    return fit_res


def fit_NB(spark, normal_matrix, index_name, output_dir, path_normal_matrix_segm, libsize_n, cores):  #TODO decide on statistical mode future
    ''' Fits negative binomial on kmers expression with DESeq2
    :param spark: spark context
    :param normal_matrix: dataframe normal matrix
    :param index_name: str kmer column name
    :param output_dir: str output directory
    :param path_normal_matrix_segm: str path to save matrix
    :param libsize_n: libsize matrix
    :param cores: int number of cores to run DESeq2
    :return: spark session, dataframe normal matrix with additional fit columns
    '''

    def nb_cdf(mean_, disp):
        probA = disp / (disp + mean_)
        N = (probA * mean_) / (1 - probA)
        return float(scipy.stats.nbinom.cdf(0, n=N, p=probA))

    design_matrix = pd.DataFrame([1] * (len(normal_matrix.columns) - 1), columns=["design"])
    design_matrix['sample'] = [col_name for col_name in normal_matrix.schema.names if col_name != index_name]

    # Run DESEq2
    design_matrix = design_matrix.set_index('sample')
    logging.info("Run DESeq2")
    normal_matrix = DESeq2(normal_matrix.toPandas().set_index(index_name), design_matrix, normalize=libsize_n, cores=cores)
    save_pd_toparquet(os.path.join(output_dir, os.path.basename(path_normal_matrix_segm).split('.')[0] + 'deseq_fit' + '.pq.gz'),
                      normal_matrix, compression = 'gzip', verbose = False)
    # Test on probability of noise
    logging.info("Test if noise")
    normal_matrix = spark.createDataFrame(normal_matrix.reset_index())
    stat_udf = sf.udf(nb_cdf, st.FloatType())
    logging.info("Filter out noise from normals")
    normal_matrix = normal_matrix.withColumn("proba_zero", stat_udf(sf.col('baseMean'), sf.col('dispersion')))

    return spark, normal_matrix


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


def pq_WithRenamedCols(spark, list_paths):
    '''
    Load parquet file, convert column names and convert to spark dataframe
    :param spark: spark context
    :param list_paths: list of paths to be read as one parquet file
    :return: loaded spark dataframe
    '''
    df = spark.read.parquet(*list_paths , mergeSchema=True)
    old_name = df.columns
    new_names = [name_.replace('-', '').replace('.', '').replace('_', '')  for name_ in  old_name]
    df = df.toDF(*new_names)
    return df


def process_matrix_file(spark, index_name, jct_col, jct_annot_col, rf_annot_col,
                        path_matrix, whitelist, cross_junction, annot_flag, tot_batches=None, batch_id=None):
    '''
    Preprocess samples if expression is stored in a multi-sample format [kmers] x [samples, metadata]
    Various preprocessing steps including filtering on junction status, on reading frame annotated status,
    select whitelist samples and make kmer unique
    :param spark: spark context
    :param index_name: str kmer column name
    :param jct_col: str junction column name
    :param jct_annot_col: str junction annotated flag column name
    :param rf_annot_col: str reading frame annotated flag column name
    :param path_matrix: str path for count matrix
    :param whitelist: list whitelist for samples
    :param cross_junction: bool whether the count matrix contains junction counts or segment counts
    :param annot_flag: list with the intruction codes on how to treat the reading frame and junction annotated flags
    :param tot_batches: int batch mode only, total number of batches
    :param batch_id: int batch mode only, id of the batch
    :return: Preprocessed spark dataframe
    '''

    def cast_type_dbl(matrix, name_list, index_name):
        return matrix.select(
            [sf.col(name_).cast(st.DoubleType()).alias(name_) if name_ != index_name else sf.col(name_) for name_ in
             name_list])

    if path_matrix is not None:
        # Rename columns for better spark processing (spark does not support '_')
        rename = True # For development
        logging.info(f'Load {path_matrix}')
        if rename:
            logging.info("Rename")
            matrix = pq_WithRenamedCols(spark, path_matrix)
        else:
            matrix = spark.read.parquet(*path_matrix , mergeSchema=True)

        # In batch mode: Dev remove kmers for the matrix based on the hash function
        if tot_batches:
            logging.info(f'Filter foreground and background based on Hash ; making {tot_batches} batches, select batch number {batch_id}')

            matrix = matrix.filter(sf.abs(sf.hash('kmer') % tot_batches) == batch_id)


        # Filter on junction status depending on the content on the matrix
        if cross_junction:
            matrix = matrix.filter(f'{jct_col} == True')
        else:
            matrix = matrix.filter(f'{jct_col} == False')
        matrix = matrix.drop(jct_col)

        # Filter according to annotation flag
        matrix = filter_on_junction_kmer_annotated_flag(matrix, jct_annot_col, rf_annot_col, annot_flag)


        # Cast type and fill nans + Reduce samples (columns) to whitelist
        logging.info("Cast types")
        if whitelist is not None:
            whitelist = pd.read_csv(whitelist, sep='\t', header=None)[0].to_list()
            whitelist = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in whitelist]
            whitelist.append(index_name)
            matrix = cast_type_dbl(matrix, whitelist, index_name)
        else:
            matrix = cast_type_dbl(matrix, matrix.schema.names, index_name)

        # Fill Nans
        logging.info("Remove Nans")
        matrix = matrix.na.fill(0)

        # Remove kmers present in the table but absents from all samples
        logging.info("Remove non expressed kmers SQL-")
        partitions_ = matrix.rdd.getNumPartitions()
        logging.info(f'...partitions: {partitions_}')
        not_null = ' OR '.join(
            [f'({col_name} != 0.0)'
             for col_name in matrix.schema.names
             if col_name not in [index_name, jct_annot_col, rf_annot_col]])  # SQL style  # All zeros
        matrix = matrix.filter(not_null)

        # Make unique on kmers
        logging.info("Make unique")
        partitions_ = matrix.rdd.getNumPartitions()
        logging.info(f'...partitions: {partitions_}')
        exprs = [sf.max(sf.col(name_)).alias(name_) for name_ in matrix.schema.names if name_ != index_name]
        matrix = matrix.groupBy(index_name).agg(*exprs) #Max operation also performed on the annot flags
        partitions_ = matrix.rdd.getNumPartitions()
        logging.info(f'...partitions: {partitions_}')
        return matrix
    else:
        return None


def filter_on_junction_kmer_annotated_flag(matrix, jct_annot_col, rf_annot_col, annot_flag):
    '''
    Filters according to the junction and kmer annotated flag
    :param matrix: spark dataframe matrix to filter
    :param jct_annot_col: string junction is annotated column
    :param rf_annot_col: string reading frame is annotated column
    :param annot_flag: list of codes for filtering 0: do nothing,
    1: keep junction annotated AND reading frame annotated,
    2: keep junction annotated AND reading frame not annotated,
    3: keep junction not annotated AND reading frame annotated,
    4: keep junction not annotated AND reading frame not annotated.
    Several codes can be provided as input  and will be combined with OR operation.
    :return: spark dataframe expression matrix filtered on junction annotated and reading frame annotated flag.
    '''
    # Keep k-mers according to annotation flag
    code_to_flag = {1: (1, 1), 2: (1, 0), 3: (0, 1), 4: (0, 0)}  # junction and reading frame annotated respectively
    flag_condition = [f'({jct_annot_col} == {code_to_flag[code][0]} AND {rf_annot_col} == {code_to_flag[code][1]})'
                      for code in annot_flag if code]
    if flag_condition:
        matrix = matrix.filter(' OR '.join(flag_condition))
    return matrix


def combine_normals(normal_segm, normal_junc, index_name):
    ''' Given two matrices [kmers] x [samples, metadata] containing segment expression and junction expression,
    takes the maximal expression for each kmer across the two matrices
    :param normal_segm: spark dataframe matrix with segment expression counts
    :param normal_junc: spark dataframe matrix with junction expression counts
    :param index_name: str name of the kmer column
    :return: spark dataframe expression matrix after combining junction and segment expression
    '''
    if (normal_segm is not None) and (normal_junc is not None):
        logging.info("Combine segment and edges")
        normal_matrix = normal_segm.union(normal_junc)
        # Take max expression between edge or segment expression
        exprs = [sf.max(sf.col(name_)).alias(name_) for name_ in normal_matrix.schema.names if name_ != index_name]
        normal_matrix = normal_matrix.groupBy(index_name).agg(*exprs)
        return normal_matrix
    elif normal_segm is None:
        return normal_junc
    elif normal_junc is None:
        return normal_segm


def outlier_filtering(normal_matrix, index_name, libsize_n, expr_high_limit_normal):
    '''
    Remove very highly expressed kmers / expression outliers before fitting DESeq2. These kmers do not follow a NB,
    besides no hypothesis testing is required to set their expression status to True
    :param normal_matrix: spark dataframe for normal matrix
    :param index_name: str kmer column name
    :param libsize_n: dataframe for library size
    :param expr_high_limit_normal: float threshold value on number of normalized reads
    :return: Two different filtered instances of spark dataframe
    '''

    # With libsize
    if libsize_n is not None:
        highly_expressed_normals = ' AND '.join(
            [f'({col_name} > { expr_high_limit_normal * libsize_n.loc[col_name, "libsize_75percent"]})'
             for col_name in normal_matrix.schema.names if col_name != index_name])  # SQL style  # Expressed kmers

        ambigous_expression_normals = ' OR '.join(
            [f'({col_name} <= { expr_high_limit_normal * libsize_n.loc[col_name, "libsize_75percent"]})'
             for col_name in normal_matrix.schema.names if col_name != index_name])  # SQL style
    # Without libsize
    else:
        highly_expressed_normals = ' AND '.join([f'({col_name} >= {expr_high_limit_normal})'
                                                 for col_name in normal_matrix.schema.names if
                                                 col_name != index_name])  # SQL style  # Expressed kmers

        ambigous_expression_normals = ' OR '.join([f'({col_name} < {expr_high_limit_normal})'
                                                   for col_name in normal_matrix.schema.names if
                                                   col_name != index_name])  # SQL style

    high_expr_normals = normal_matrix.filter(highly_expressed_normals).select(sf.col(index_name))
    normal_matrix = normal_matrix.filter(ambigous_expression_normals)  # TODO add condition empty matrix
    return high_expr_normals, normal_matrix


def filter_statistical(spark, tissue_grp_files, normal_matrix, index_name, path_normal_matrix_segm, libsize_n,  #TODO decide on statistical mode future
                       threshold_noise, output_dir, cores):
    '''
    :param spark: spark context
    :param tissue_grp_files: list paths of group tissues
    :param normal_matrix: spark dataframe for normal matrix
    :param index_name: str kmer column name
    :param path_normal_matrix_segm: str path for normal matrix with segment counts
    :param libsize_n: dataframe for library size
    :param threshold_noise: float threshold for probability of being noise (< threshold)
    :param output_dir: str output directory
    :param cores: int number of cores
    :return: nothing. Function not finished
    '''
    if tissue_grp_files is not None:
        modelling_grps = []
        for tissue_grp in tissue_grp_files:
            grp = pd.read_csv(tissue_grp, header=None)[0].to_list()
            grp = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in grp]
            grp.append(index_name)
            modelling_grps.append(grp)
        else:
            modelling_grps = [[name_ for name_ in normal_matrix.schema.names if name_ != index_name]]

        logging.info(">>>... Fit Negative Binomial distribution on normal kmers ")
        for grp in modelling_grps:
            # Fit NB and Perform hypothesis testing
            normal_matrix = fit_NB(spark, normal_matrix, index_name, output_dir, path_normal_matrix_segm,
                                   libsize_n, cores)
            normal_matrix = normal_matrix.filter(sf.col("proba_zero") < threshold_noise)  # Expressed kmers

        # Join on the kmers segments. Take the kmer which junction expression is not zero everywhere


def filter_hard_threshold(matrix, index_name, jct_annot_col, rf_annot_col, libsize,
                          out_dir, expr_limit, n_samples_lim, target_sample='', annot_flag=[],
                          tag='normals', batch_tag=''):
    '''
    Filter samples based on >0 and >=X reads expressed (expansive operations) and save intermediate files.
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
    :param out_dir: str path for output directory
    :param expr_limit: float expression limit threshold to keep a kmer
    :param n_samples_lim: int number of samples that need to pass the expression limit
    :param target_sample: str name of the sample of interest.
    To be excluded in the number of samples that pass the expression limit
    :param annot_flag: list with the intruction codes on how to treat the reading frame and junction annotated flags
    :param tag: str tag related to the type of samples. Example cancer or normal
    :param batch_tag: str batch mode, batch tag to be appened to intermediate file
    :return: path_e, path_s respectively path where
    the expression-filtered matrix and the sample-filtered matrix are saved
    '''

    path_e = None
    path_s = None
    base_n_samples = 1
    base_expr = 0.0
    if target_sample:
        suffix = f'Except{target_sample}'
    else:
        suffix = ''


    if target_sample:
        logging.info(f'Target sample {target_sample} not included in the cohort filtering')
    if libsize is not None:
        matrix = matrix.select(index_name, *[
            sf.round(sf.col(name_) / libsize.loc[name_, "libsize_75percent"], 2).alias(name_)
            for name_ in matrix.schema.names if name_ not in [index_name, jct_annot_col, rf_annot_col]])

    # Expression filtering, take k-mers with >= X reads in >= 1 sample
    if expr_limit:
        path_e = os.path.join(out_dir,f'interm_{tag}_combiExprCohortLim{expr_limit}Across{base_n_samples}{suffix}{batch_tag}' + '.tsv')
        if not os.path.isfile(os.path.join(path_e, '_SUCCESS')):
            logging.info(f'Filter matrix with cohort expression support >= {expr_limit} in {base_n_samples} sample')
            matrix_e = matrix.select(index_name, *[
                sf.when(sf.col(name_) >= expr_limit, 1).otherwise(0).alias(name_)
                for name_ in matrix.schema.names if name_
                not in [target_sample, index_name, jct_annot_col, rf_annot_col]])
            matrix_e = matrix_e.rdd.map(tuple).map(lambda x: (x[0], sum(x[1:]))
                                                                 ).filter(lambda x: x[1] >= base_n_samples)
            logging.info(f'Save to {path_e}')
            matrix_e.map(lambda x: "%s\t%s" % (x[0], x[1])).saveAsTextFile(path_e)
        else:
            logging.info(
            f'Filter matrix with cohort expression support {expr_limit} in {base_n_samples} sample(s) already performed. Loading results from {path_e}')
            logging.info(f'Using intermediate files means ignoring --annotated-flags {annot_flag} parameter.')

    # Sample filtering, take k-mers with exclude >0 reads in >= 1 sample
    if n_samples_lim is not None:
        path_s = os.path.join(out_dir, f'interm_{tag}_combiExprCohortLim{base_expr}Across{base_n_samples}{suffix}{batch_tag}' + '.tsv')
        if not os.path.isfile(os.path.join(path_s, '_SUCCESS')):
            logging.info(f'Filter matrix with cohort expression support > {base_expr} in {base_n_samples} sample')
            matrix_s = matrix.select(index_name, *[
                sf.when(sf.col(name_) > base_expr, 1).otherwise(0).alias(name_)
                for name_ in matrix.schema.names if name_
                not in [target_sample, index_name, jct_annot_col, rf_annot_col]])
            matrix_s = matrix_s.rdd.map(tuple).map(lambda x: (x[0], sum(x[1:]))
                                                                 ).filter(lambda x: x[1] >= base_n_samples)
            logging.info(f'Save to {path_s}')
            matrix_s.map(lambda x: "%s\t%s" % (x[0], x[1])).saveAsTextFile(path_s)
        else:
            logging.info(
            f'Filter matrix with cohort expression support {base_expr} in {base_n_samples} sample(s) already performed. Loading results from {path_s}')
            logging.info(f'Using intermediate files means ignoring --annotated-flags {annot_flag} parameter.')
    return path_e, path_s


def combine_hard_threshold_normals(spark, path_normal_kmers_e, path_normal_kmers_s, n_samples_lim_normal, index_name):
    '''
    Filter samples based on X reads and H samples.
    a. kmer need to be if >= X reads in >= 1 sample -> load path_*_e (expression) file
    b. kmer needs to be >0 reads in >= H samples -> load path_*_s (sample) and apply H threshold on the fly
    The kmers selected need to pass a. OR b.
    :param spark: spark context
    :param path_normal_kmers_e: str path for intermediate file which applied >= X reads in >= 1 sample
    :param path_normal_kmers_s: str path for intermediate file which applied >0 reads in >= 1 sample
    :param n_samples_lim_normal: int number of samples in which any number of reads is required
    :param index_name: index_name: str kmer column name
    :return: spark dataframe filtered with >= X reads in 1 sample OR >0 reads in H samples
    '''
    if path_normal_kmers_e: #a.k.a exclude >= X reads in >= 1 sample)
        normal_matrix_e = spark.read.csv(path_normal_kmers_e, sep=r'\t', header=False)
        normal_matrix_e = normal_matrix_e.withColumnRenamed('_c0', index_name)
        normal_matrix_e = normal_matrix_e.select(sf.col(index_name))
        if not path_normal_kmers_s:
            return normal_matrix_e
    if path_normal_kmers_s:  # (a.k.a exclude >0  reads in >= H samples)
        normal_matrix_s = spark.read.csv(path_normal_kmers_s, sep=r'\t', header=False)
        normal_matrix_s = normal_matrix_s.withColumnRenamed('_c0', index_name)
        normal_matrix_s = normal_matrix_s.withColumnRenamed('_c1', "n_samples")
        if n_samples_lim_normal > 1:
            logging.info( f'Filter matrix with cohort expression support > {0} in {n_samples_lim_normal} sample(s)')
            normal_matrix_s = normal_matrix_s.filter(sf.col('n_samples') >= n_samples_lim_normal)
        normal_matrix_s = normal_matrix_s.select(sf.col(index_name))
        if not path_normal_kmers_e:
            return normal_matrix_s
    if path_normal_kmers_e and path_normal_kmers_s:
        normal_matrix_res = normal_matrix_e.union(normal_matrix_s) # Do not make distinct
        return normal_matrix_res
    else:
        return None


def combine_hard_threshold_cancers(spark, cancer_matrix, path_cancer_kmers_e, path_cancer_kmers_s,
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
    :param cohort_expr_support_cancer: float expression threshold to be met for the kmer
    :param n_samples_lim_cancer: int number of samples in which the expression threshold need to be met
    :param index_name: str kmer column name
    :return: spark dataframe with foreground filtered for >= X reads and in >= H samples
    '''
    # Load intermediate file
    if cohort_expr_support_cancer != 0.0:
        valid_foreground = spark.read.csv(path_cancer_kmers_e, sep=r'\t', header=False)
    else:
        valid_foreground = spark.read.csv(path_cancer_kmers_s, sep=r'\t', header=False)
    valid_foreground = valid_foreground.withColumnRenamed('_c0', index_name)
    valid_foreground = valid_foreground.withColumnRenamed('_c1', "n_samples")
    # kmer need to be expressed with >= X reads and this in >= H samples
    if n_samples_lim_cancer > 1:
        logging.info( f'Filter matrix with cohort expression support >= {cohort_expr_support_cancer} in {n_samples_lim_cancer} sample(s)')
        valid_foreground = valid_foreground.filter(sf.col('n_samples') >= n_samples_lim_cancer)
    valid_foreground = valid_foreground.select(sf.col(index_name))
    # Apply to preprocessed matrix
    cancer_cross_filter = cancer_matrix.join(valid_foreground, ["kmer"],
                                             how='right')
    return cancer_cross_filter


def preprocess_kmer_file(cancer_kmers, cancer_sample, drop_cols, expression_fields,
                         jct_col, jct_annot_col, rf_annot_col, index_name, libsize_c, annot_flag, cross_junction):
    '''
    Preprocess samples if expression is stored in a single sample file:
    [kmers] x [junction expression, segment expression, metadata]
    Performs various preprocessing steps including filtering on junction status and on reading frame annotated status,
    make kmer unique, normalize
    :param cancer_kmers: spark dataframe with expression counts for cancer
    :param cancer_sample: str cancer sample ID
    :param drop_cols: list columns to be dropped
    :param expression_fields: list containing expression column names
    :param jct_col: str junction column name
    :param jct_annot_col: str junction annotated flag column name
    :param rf_annot_col: str reading frame annotated flag column name
    :param index_name: str kmer column name
    :param libsize_c: dataframe with library size
    :param annot_flag: list with the intruction codes on how to treat the reading frame and junction annotated flags
    :param cross_junction bool whether the count matrix contains junction counts or segment counts
    :return: spark dataframe for preprocessed cancer kmers matrix
    '''

    def collapse_values(value):
        return max([np.float(i) if i != 'nan' else 0.0 for i in value.split('/')])  # np.nanmax not supported

    # Filter on juction status
    if cross_junction == 1:
        cancer_kmers = cancer_kmers.filter(f'{jct_col} == True')
    elif cross_junction == 0:
        cancer_kmers = cancer_kmers.filter(f'{jct_col} == False')

    # Drop junction column
    for drop_col in drop_cols:
        cancer_kmers = cancer_kmers.drop(sf.col(drop_col))

    # Keep k-mers according to annotation flag
    cancer_kmers = filter_on_junction_kmer_annotated_flag(cancer_kmers, jct_annot_col, rf_annot_col, annot_flag)

    logging.info("Collapse kmer horizontal")
    # Remove the '/' in the expression data (kmers duplicate within a gene have 'expression1/expression2' format
    local_max = sf.udf(collapse_values, st.FloatType())
    for name_ in expression_fields:
        cancer_kmers = cancer_kmers.withColumn(name_, local_max(name_))

    # Make kmers unique (Take max expression)
    logging.info("Collapse kmer vertical")
    cancer_kmers = cancer_kmers.withColumn(jct_col, sf.col(jct_col).cast("boolean").cast("int"))
    exprs = [sf.max(sf.col(name_)).alias(name_) for name_ in cancer_kmers.schema.names if name_ != index_name]
    cancer_kmers = cancer_kmers.groupBy(index_name).agg(*exprs)

    # Remove kmers unexpressed (both junction and segment expression null)
    cancer_kmers = cancer_kmers.withColumn('allnull', sum(cancer_kmers[name_] for name_ in expression_fields))
    cancer_kmers = cancer_kmers.filter(sf.col("allnull") > 0.0)
    cancer_kmers = cancer_kmers.drop("allnull")

    # Normalize by library size
    if libsize_c is not None:
        for name_ in expression_fields:
            cancer_kmers = cancer_kmers.withColumn(name_, sf.round(
                cancer_kmers[name_] / libsize_c.loc[cancer_sample, "libsize_75percent"], 2))
    else:
        for name_ in expression_fields:
            cancer_kmers = cancer_kmers.withColumn(name_, sf.round(cancer_kmers[name_], 2))

    return cancer_kmers


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
        # if  max( edge expression 1 and 2)<threshold and  max( segment expression 1 and 2)>= threshold: keep
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


def save_spark(cancer_kmers, output_dir, path_final_fil, outpartitions=None):
    '''
    Saves a spark datadrame as a single or partitionned csv file
    :param cancer_kmers: spark dataframe matrix with expression counts for cancer
    :param output_dir: str path for output directory
    :param path_final_fil: str path to save the spark dataframe
    :param outpartitions: int number of partitions for saving
    '''
    # save
    logging.info(f'Save to {path_final_fil}')
    pathlib.Path(output_dir).mkdir(exist_ok=True, parents=True)
    if outpartitions is not None:
        cancer_kmers.repartition(outpartitions).write.mode('overwrite').options(header="true", sep="\t").csv(path_final_fil)
    else:
        cancer_kmers.write.mode('overwrite').options(header="true", sep="\t").csv(path_final_fil)


def loader(spark, path_kmer_list):
    '''
    Loads path parquet or csv to spark dataframe
    If user provides a list of parquet files with the same schema -> They will be loaded as a single file
    elif user provides a list with a single csv -> Only the first entry of the list will be loaded as a single file
    :param spark: spark context
    :param path_kmer_list: list of paths to read
    :return: loaded spark dataframe
    '''
    #TODO allow multiple tsv
    if 'tsv' in path_kmer_list[0]:
        logging.warning(f'Only the first file of {path_kmer_list} will be read. Use list of parquets to process multiple paths')
        matrix = spark.read.csv(path_kmer_list[0], sep=r'\t', header=False)
    else:
        matrix = spark.read.parquet(*path_kmer_list)

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
        header = f'{"sample"}\t{"mutation_mode"}\t{"min_sample_reads"}\t{"#_of_cohort_samples"}\t{"reads_per_cohort_sample"}\t{"#_normal_samples_allowed"}\t{"normal_cohort_id"}\t{"reads_per_normal_sample"}'
        line =   f'{cancer_sample_ori}\t{mutation_mode}\t{sample_expr_support_cancer}\t{n_samples_lim_cancer}\t{cohort_expr_support_cancer}\t{n_samples_lim_normal}\t{id_normals}\t{cohort_expr_support_normal}'
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


def redirect_scratch(scratch_dir, interm_dir_norm, interm_dir_canc, output_dir):
    '''
    Set the directory to save intermediary file
    - Uses either the scratch directory variable from the cluster
    - The output directory
    - Any other specified normal or cancer scratch directory
    Default. Uses output directory
    :param scratch_dir: str os environement variable name containing the cluster scratch directory path
    :param interm_dir_norm: str custom scatch dir path to save intermediate normal files
    :param interm_dir_canc: str custom scatch dir path to save intermediate cancer files
    :param output_dir: str output directory for the filtered matrix
    :return:
    '''
    if scratch_dir:
        cancer_out = os.environ[scratch_dir]
        normal_out = os.environ[scratch_dir]
        return normal_out, cancer_out
    if interm_dir_canc:
        cancer_out = interm_dir_canc
    else:
        cancer_out = output_dir
    if  interm_dir_norm:
        normal_out = interm_dir_norm
    else:
        normal_out = output_dir
    return normal_out, cancer_out
