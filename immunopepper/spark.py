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


def process_libsize(path_lib, custom_normalizer):
    lib = pd.read_csv(path_lib, sep='\t')
    if custom_normalizer:
        lib['libsize_75percent'] = lib['libsize_75percent'] / custom_normalizer
    else:
        lib['libsize_75percent'] = lib['libsize_75percent'] / np.median(lib['libsize_75percent'])
    lib['sample'] = [sample.replace('-', '').replace('.', '').replace('_','') for sample in lib['sample']]
    lib = lib.set_index('sample')
    return lib

def pq_WithRenamedCols(spark, list_paths, outdir):
    df = spark.read.parquet(*list_paths , mergeSchema=True)
    old_name = df.columns
    new_names = [name_.replace('-', '').replace('.', '').replace('_', '')  for name_ in  old_name]
    df = df.toDF(*new_names)
    return df


def process_matrix_file(spark, index_name, jct_col, path_normal_matrix, outdir, whitelist, parallelism, cross_junction, tot_batches=None, batch_id=None):
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
        rename = True # For development
        logging.info("Load {}".format(path_normal_matrix))
        if rename:
            logging.info("Rename")
            normal_matrix = pq_WithRenamedCols(spark, path_normal_matrix, outdir)
        else:
            normal_matrix = spark.read.parquet(*path_normal_matrix , mergeSchema=True)

        # Dev remove kmers for the matrix based on the hash function
        if tot_batches:
            logging.info("Filter foreground and background based on Hash ; making {} batches, select batch number {}".format(tot_batches, batch_id))

            normal_matrix = normal_matrix.filter(sf.abs(sf.hash('kmer') % tot_batches) == batch_id)


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

def filter_hard_threshold(normal_matrix, index_name, libsize, out_dir, expr_limit, n_samples_lim, target_sample='', tag='normals', batch_tag=''):
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
    path_e = None
    path_s = None
    base_n_samples  = 1
    base_expr = 0.0 
    if target_sample:
        suffix = 'Except{}'.format(target_sample)
    else:
        suffix = ''


    if target_sample:
        logging.info("Target sample {} not included in the cohort filtering".format(target_sample))
    if libsize is not None:
        normal_matrix = normal_matrix.select(index_name, *[
            sf.round(sf.col(name_) / libsize.loc[name_, "libsize_75percent"], 2).alias(name_)
            for name_ in normal_matrix.schema.names if name_ != index_name])

    if expr_limit:  #a.k.a exclude >= X reads in >= 1 sample)
        path_e = os.path.join(out_dir,'interm_{}_combiExprCohortLim{}Across{}{}{}'.format(tag, expr_limit, base_n_samples, suffix, batch_tag) + '.tsv')
        if not os.path.isfile(os.path.join(path_e, '_SUCCESS')):
            logging.info("Filter matrix with cohort expression support >= {} in {} sample".format(expr_limit, base_n_samples))
            normal_matrix_e = normal_matrix.select(index_name, *[
                sf.when(sf.col(name_) >= expr_limit, 1).otherwise(0).alias(name_)
                for name_ in normal_matrix.schema.names if (name_ != index_name) and (name_ != target_sample ) ])
            normal_matrix_e = normal_matrix_e.rdd.map(tuple).map(lambda x: (x[0], sum(x[1:]))).filter(lambda x: x[1] >= base_n_samples)
            logging.info("Save to {}".format(path_e))
            normal_matrix_e.map(lambda x: "%s\t%s" % (x[0], x[1])).saveAsTextFile(path_e)
        else:
            logging.info(
            "Filter matrix with cohort expression support {} in {} sample(s) already performed. Loading results from {}".format(
                expr_limit, base_n_samples, path_e))
            
    if n_samples_lim is not None: # (a.k.a exclude >0  reads in >= H samples) --> H samples filtering done subsequently # n_samples_lim can be 0 -> 1 i used
        path_s = os.path.join(out_dir,'interm_{}_combiExprCohortLim{}Across{}{}{}'.format(tag, base_expr, base_n_samples, suffix, batch_tag) + '.tsv')
        if not os.path.isfile(os.path.join(path_s, '_SUCCESS')):
            logging.info("Filter matrix with cohort expression support > {} in {} sample".format(base_expr, base_n_samples))
            normal_matrix_s = normal_matrix.select(index_name, *[
                sf.when(sf.col(name_) > base_expr, 1).otherwise(0).alias(name_)
                for name_ in normal_matrix.schema.names if (name_ != index_name) and (name_ != target_sample ) ])
            normal_matrix_s = normal_matrix_s.rdd.map(tuple).map(lambda x: (x[0], sum(x[1:]))).filter(lambda x: x[1] >= base_n_samples)
            logging.info("Save to {}".format(path_s))
            normal_matrix_s.map(lambda x: "%s\t%s" % (x[0], x[1])).saveAsTextFile(path_s)
        else:
            logging.info(
            "Filter matrix with cohort expression support {} in {} sample(s) already performed. Loading results from {}".format(
                base_expr, base_n_samples, path_s))

    return path_e, path_s


def combine_hard_threshold_normals(spark, path_normal_kmers_e, path_normal_kmers_s, n_samples_lim_normal, index_name):
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
            logging.info( "Filter matrix with cohort expression support > {} in {} sample(s)".format(0, n_samples_lim_normal))
            normal_matrix_s = normal_matrix_s.filter(sf.col('n_samples') >= n_samples_lim_normal)
        normal_matrix_s = normal_matrix_s.select(sf.col(index_name))
        if not path_normal_kmers_e:
            return normal_matrix_s
    if path_normal_kmers_e and path_normal_kmers_s:
        normal_matrix_res = normal_matrix_e.union(normal_matrix_s) # Do not make distinct
        return normal_matrix_res
    else:
        return None


def combine_hard_threshold_cancers(spark, cancer_matrix, path_cancer_kmers_e, cohort_expr_support_cancer, n_samples_lim_cancer, index_name, cancer_sample):
    valid_foreground = spark.read.csv(path_cancer_kmers_e, sep=r'\t', header=False)
    valid_foreground = valid_foreground.withColumnRenamed('_c0', index_name)
    valid_foreground = valid_foreground.withColumnRenamed('_c1', "n_samples")
    if n_samples_lim_cancer > 1:
        logging.info( "Filter matrix with cohort expression support >= {} in {} sample(s)".format(cohort_expr_support_cancer, n_samples_lim_cancer))
        valid_foreground = valid_foreground.filter(sf.col('n_samples') >= n_samples_lim_cancer)
    valid_foreground = valid_foreground.select(sf.col(index_name))
    cancer_cross_filter = cancer_matrix.join(valid_foreground, ["kmer"],  # Probably do union differently
                                             how='right').select([index_name, cancer_sample]) # Intersect with preprocessed cancer matrix
    return cancer_cross_filter

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

    logging.info("...partitions: {}".format(matrix_kmers.rdd.getNumPartitions()))
    if threshold != 0:
        logging.info("Filter with {} >= {}".format(filter_field, threshold))
        matrix_kmers = matrix_kmers.filter(sf.col(filter_field) >= threshold)
    else:
        logging.info("Filter with {} > {}".format(filter_field, threshold))
        matrix_kmers = matrix_kmers.filter(sf.col(filter_field) > threshold)
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



def loader(spark, path_normal_kmer_list):
    #TODO allow multiple tsv
    if 'tsv' in path_normal_kmer_list[0]:
        logging.warning("Only the first file of {} will be read. Use list of parquets to process multiple paths".format(path_normal_kmer_list))
        normal_matrix = spark.read.csv(path_normal_kmer_list[0], sep=r'\t', header=False)
    else:
        normal_matrix = spark.read.parquet(*path_normal_kmer_list)

    return normal_matrix


def output_count(perform_count, matrix, report_count, report_step, step_string):
    if perform_count:
        mycount = matrix.count()
        report_count.append(mycount)
        report_step.append(step_string)
        logging.info('# {} n = {} kmers'.format(step_string, mycount))

def save_output_count(output_count, report_count, report_steps, prefix, cancer_sample_ori, mutation_mode,
                      sample_expr_support_cancer, cohort_expr_support_cancer, n_samples_lim_cancer,
                          cohort_expr_support_normal, n_samples_lim_normal, id_normals):
    if output_count:
        header = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format("sample", "mutation_mode", "min_sample_reads", "#_of_cohort_samples", "reads_per_cohort_sample", "#_normal_samples_allowed", "normal_cohort_id", "reads_per_normal_sample")
        line =   '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(cancer_sample_ori, mutation_mode, sample_expr_support_cancer, n_samples_lim_cancer, cohort_expr_support_cancer, n_samples_lim_normal, id_normals, cohort_expr_support_normal)
        for idx in np.arange(len(report_count)):
            header += "\t{}".format(report_steps[idx])
            line += "\t{}".format(report_count[idx])
        if prefix:
            header += "\t{}".format("info")
            line += "\t{}".format(prefix)
        header += "\n"
        line += "\n"
        if not os.path.exists(output_count):
            with open(output_count,"w") as f:
                f.write(header)
        with open(output_count, "a") as f:
            f.write(line)
        logging.info("Save intermediate info to {}".format(output_count))

def redirect_scratch(scratch_dir, interm_dir_norm, interm_dir_canc, output_dir):
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
