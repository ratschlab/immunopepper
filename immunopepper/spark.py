import logging
import numpy as np
import os
import pandas as pd
import pathlib
import pyarrow.parquet as pq
from pyspark.sql import functions as sf
from pyspark.sql import types as st
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, Formula, r
from rpy2.robjects.packages import importr
import scipy
from scipy import stats
import sys

from .io_ import save_pd_toparquet

pandas2ri.activate()
deseq = importr('DESeq2')
BiocParallel = importr('BiocParallel')
BiocGenerics = importr("BiocGenerics")



def DESeq2(count_matrix, design_matrix, normalize, cores=1):
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
        dds0.do_slot('colData').do_slot('listData')[1] = ro.vectors.FloatVector(list(normalize.loc[order_size_factor, 'libsize_75percent'])) # Enforce size factors
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


def fit_NB(spark, normal_matrix, index_name, output_dir, path_normal_matrix_segm, libsize_n, cores):
    ''' Fits negative binomial on kmers expression with DESeq2
    Parameters
    ---------
    spark: spark session
    normal_matrix: normal matrix
    index_name: kmer column name
    output_dir: output directory
    path_normal_matrix_segm: path to save matrix
    libsize_n: libsize matrix
    cores: number of cores to run DESeq2
    Returns
    --------
    spark: spark session
    normal_matrix: normal matrix
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


def process_libsize(path_lib):
    lib = pd.read_csv(path_lib, sep='\t')
    lib['libsize_75percent'] = lib['libsize_75percent'] / np.median(lib['libsize_75percent'])
    lib['sample'] = [sample.replace('-', '').replace('.', '').replace('_','') for sample in lib['sample']]
    lib = lib.set_index('sample')
    return lib

def pq_WithRenamedCols(list_paths, outdir):
    list_path_tmp = []
    for path in list_paths:
        df = pq.read_table(path)
        df = df.rename_columns(name_.replace('-', '').replace('.', '').replace('_', '') for name_ in df.schema.names) # characters causing issue in spark
        path_tmp = os.path.join(outdir, os.path.basename(path).split('.')[0] + '_tmp' + '.pq')
        if os.path.exists(path_tmp):
            os.remove(path_tmp)
        pqwriter = pq.ParquetWriter(path_tmp, df.schema, compression=None)
        pqwriter.write_table(df)
        pqwriter.close()
        list_path_tmp.append(path_tmp)
    return list_path_tmp


def process_matrix_file(spark, index_name, jct_col, path_normal_matrix, outdir, whitelist, parallelism, cross_junction):
    ''' Preprocess normal samples
    - corrects names
    - corrects types
    - make unique
    Parameters:
    ----------
    spark: spark context
    index_name: kmer column name
    jct_col: junction column name
    path_normal_matrix_segm: path for normal matrix
    whitelist: whitelist for normal samples
    Returns :
    ----------
    spark: spark session
    normal_matrix: Preprocessed normal matrix
        '''

    def cast_type_dbl(normal_matrix, name_list, index_name):
        return normal_matrix.select(
            [sf.col(name_).cast(st.DoubleType()).alias(name_) if name_ != index_name else sf.col(name_) for name_ in
             name_list])

    if path_normal_matrix is not None:
        # Rename
        rename = True  # For development
        logging.info("Load {}".format(path_normal_matrix))
        if rename:
            logging.info("Rename")
            path_normal_matrix_tmp = pq_WithRenamedCols(path_normal_matrix, outdir)
            normal_matrix = spark.read.parquet(*path_normal_matrix_tmp)
        else:
            normal_matrix = spark.read.parquet(*path_normal_matrix)

        # Keep relevant junction status and drop junction column
        if cross_junction:
            normal_matrix = normal_matrix.filter("{} == True".format(jct_col))
        else:
            normal_matrix = normal_matrix.filter("{} == False".format(jct_col))
        normal_matrix = normal_matrix.drop(jct_col)

        # Cast type and fill nans + Reduce samples (columns) to whitelist
        logging.info("Cast types")
        if whitelist is not None:
            whitelist = pd.read_csv(whitelist, sep='\t', header=None)[0].to_list()
            whitelist = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in whitelist]
            whitelist.append(index_name)
            normal_matrix = cast_type_dbl(normal_matrix, whitelist, index_name)
        else:
            normal_matrix = cast_type_dbl(normal_matrix, normal_matrix.schema.names, index_name)

        # Fill Nans
        logging.info("Remove Nans")
        normal_matrix = normal_matrix.na.fill(0)

        # Remove kmers abscent from all samples
        logging.info("Remove non expressed kmers SQL-")
        logging.info("...partitions: {}".format(normal_matrix.rdd.getNumPartitions()))
        not_null = ' OR '.join(
            ['({} != 0.0)'.format(col_name)
             for col_name in normal_matrix.schema.names if col_name != index_name])  # SQL style  # All zeros
        normal_matrix = normal_matrix.filter(not_null)

        # Make unique
        logging.info("Make unique")
        logging.info("...partitions: {}".format(normal_matrix.rdd.getNumPartitions()))
        exprs = [sf.max(sf.col(name_)).alias(name_) for name_ in normal_matrix.schema.names if name_ != index_name]
        normal_matrix = normal_matrix.groupBy(index_name).agg(*exprs)
        logging.info("...partitions: {}".format(normal_matrix.rdd.getNumPartitions()))
        return normal_matrix
    else:
        return None



def combine_normals(normal_segm, normal_junc, index_name):
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
    ''' Remove very highly expressed kmers / expression outliers before fitting DESeq2. These kmers do not follow a NB,
    besides no hypothesis testing is required to set their expression status to True
    Parameters:
    spark: spark context
    normal_matrix: normal matrix
    index_name: kmer column name
    libsize_n: libsize matrix
    expr_high_limit_normal: normalized count limit for highly expressed kmers

    Returns :
    ----------
    spark: spark session
    normal_matrix: Preprocessed normal matrix

    '''

    # With libsize
    if libsize_n is not None:
        highly_expressed_normals = ' AND '.join(
            ['({} > {})'.format(col_name, expr_high_limit_normal * libsize_n.loc[col_name, "libsize_75percent"])
             for col_name in normal_matrix.schema.names if col_name != index_name])  # SQL style  # Expressed kmers

        ambigous_expression_normals = ' OR '.join(
            ['({} <= {})'.format(col_name, expr_high_limit_normal * libsize_n.loc[col_name, "libsize_75percent"])
             for col_name in normal_matrix.schema.names if col_name != index_name])  # SQL style
    # Without libsize
    else:
        highly_expressed_normals = ' AND '.join(['({} >= {})'.format(col_name, expr_high_limit_normal)
                                                 for col_name in normal_matrix.schema.names if
                                                 col_name != index_name])  # SQL style  # Expressed kmers

        ambigous_expression_normals = ' OR '.join(['({} < {})'.format(col_name, expr_high_limit_normal)
                                                   for col_name in normal_matrix.schema.names if
                                                   col_name != index_name])  # SQL style

    high_expr_normals = normal_matrix.filter(highly_expressed_normals).select(sf.col(index_name))
    normal_matrix = normal_matrix.filter(ambigous_expression_normals)  # TODO add condition empty matrix
    return high_expr_normals, normal_matrix


def filter_statistical(spark, tissue_grp_files, normal_matrix, index_name, path_normal_matrix_segm, libsize_n,
                       threshold_noise, output_dir, cores):
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

def filter_hard_threshold(normal_matrix, index_name, libsize, out_dir, expr_limit, n_samples_lim, target_sample='', tag='normals' ):
    ''' Filter samples based on j reads in at least n samples. The expressions are normalized for library size
    The filtering is either performed on the full cohort in the matrix or in the cohort excluding a target sample

    Parameters:
    ----------
    spark session
    normal matrix
    index column name
    libsize matrix
    expr_limit (j reads)
    n_samples_lim (n samples)
    target_sample
    Returns :
    ----------
    spark context
    Filtered normal matrix
        '''

    logging.info("Filter matrix with cohort expression support {} in {} samples".format(expr_limit, n_samples_lim))
    if target_sample:
        logging.info("Target sample {} not included in the cohort filtering".format(target_sample))
    if libsize is not None:
        normal_matrix = normal_matrix.select(index_name, *[
            sf.round(sf.col(name_) / libsize.loc[name_, "libsize_75percent"], 2).alias(name_)
            for name_ in normal_matrix.schema.names if name_ != index_name])

    normal_matrix = normal_matrix.select(index_name, *[
        sf.when(sf.col(name_) >= expr_limit, 1).otherwise(0).alias(name_)
        for name_ in normal_matrix.schema.names if (name_ != index_name) and (name_ != target_sample ) ]) #TODO TEST LINE

    normal_matrix = normal_matrix.rdd.map(tuple).map(lambda x: (x[0], sum(x[1:]))).filter(lambda x: x[1] >= n_samples_lim)

    if target_sample:
        suffix = 'Except{}'.format(target_sample)
    else:
        suffix = ''
    path_ = os.path.join(out_dir,
                         'interm_{}_combiExprCohortLim{}Across{}{}'.format( tag, expr_limit,
                             n_samples_lim, suffix) + '.tsv')
    logging.info("Save to {}".format(path_))
    normal_matrix.map(lambda x: "%s\t%s" % (x[0], x[1])).saveAsTextFile(path_)

    return path_


def preprocess_kmer_file(cancer_kmers, cancer_sample, drop_cols, expression_fields, jct_col, index_name, libsize_c, cross_junction):
    ''' Preprocess cancer samples
    - Make kmers unique
    - Filter kmers on junction status
    - Normalize
    Parameters:
    ----------
    cancer_kmers: cancer kmer matrix
    cancer_sample: associated cancer ID
    drop_cols: colums to be dropped
    expression_fields: list of segment and junction expression column names
    jct_col: junction status column name
    index_name: kmer column name
    libsize_c: libsize matrix for cancer samples
    cross_junction: Information to filter on juction status. None (both, no filtering), True (junction), False (non junction)
    Returns
    --------
    cancer_kmers: cancer kmers matrix,
    cancer_path_tmp: path of renamed temporary file
    jct_type: string indicating which junction filtering has been performed
    '''

    def collapse_values(value):
        return max([np.float(i) if i != 'nan' else 0.0 for i in value.split('/')])  # np.nanmax not supported

    # Filter on juction status
    if cross_junction == 1:
        cancer_kmers = cancer_kmers.filter("{} == True".format(jct_col))
    elif cross_junction == 0:
        cancer_kmers = cancer_kmers.filter("{} == False".format(jct_col))

    # Drop junction column
    for drop_col in drop_cols:
        cancer_kmers = cancer_kmers.drop(sf.col(drop_col))
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


def filter_expr_kmer(matrix_kmers, filter_field, threshold):
    logging.info("Filter out if {} <= {}".format(filter_field, threshold))
    logging.info("...partitions: {}".format(matrix_kmers.rdd.getNumPartitions()))
    matrix_kmers = matrix_kmers.filter(sf.col(filter_field) >= threshold)
    return matrix_kmers


def combine_cancer(cancer_kmers_segm, cancer_kmers_edge, index_name):
    if (cancer_kmers_segm is not None) and (cancer_kmers_edge is not None):
        logging.info("Combine segment and edges")
        cancer_kmers_segm = cancer_kmers_segm.join(cancer_kmers_edge,
                                                   cancer_kmers_segm[index_name] == cancer_kmers_edge[index_name],
                                                   how='left_anti')   # if  max( edge expression 1 and 2)<threshold and  max( segment expression 1 and 2)>= threshold: keep
        cancer_kmers = cancer_kmers_edge.union(cancer_kmers_segm)
        logging.info("...partitions cancer filtered: {}".format(cancer_kmers.rdd.getNumPartitions()))
        return cancer_kmers
    elif (cancer_kmers_segm is None):
        return cancer_kmers_edge
    elif (cancer_kmers_edge is None):
        return cancer_kmers_segm


def remove_uniprot(spark, cancer_kmers, uniprot, index_name):
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
    # save
    logging.info("Save to {}".format(path_final_fil))
    pathlib.Path(output_dir).mkdir(exist_ok=True, parents=True)
    if outpartitions is not None:
        cancer_kmers.repartition(outpartitions).write.mode('overwrite').options(header="true",sep="\t").csv(path_final_fil)
    else:
        cancer_kmers.write.mode('overwrite').options(header="true",sep="\t").csv(path_final_fil)



