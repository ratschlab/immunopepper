from functools import reduce
import logging
import numpy as np
import os
from operator import add
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
import shutil
import sys


from .config import create_spark_session_from_config
from .config import default_spark_config
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

    return spark , normal_matrix


def process_libsize(path_lib):
    lib = pd.read_csv(path_lib, sep='\t')
    lib['libsize_75percent'] = lib['libsize_75percent'] / np.median(lib['libsize_75percent'])
    lib['sample'] = [sample.replace('-', '').replace('.', '').replace('_','') for sample in lib['sample']]
    lib = lib.set_index('sample')
    return lib

def pq_WithRenamedCols(path):
    df = pq.read_table(path)
    df = df.rename_columns(name_.replace('-', '').replace('.', '').replace('_', '') for name_ in df.schema.names) # characters causing issue in spark
    path_tmp = path + 'tmp'
    if os.path.exists(path_tmp):
        os.remove(path_tmp)
    pqwriter = pq.ParquetWriter(path_tmp, df.schema, compression=None)
    pqwriter.write_table(df)
    pqwriter.close()
    return path_tmp


def process_normals(spark, index_name, jct_col, path_normal_matrix_segm, whitelist):
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
    # Rename
    rename = False  # For development
    if rename:
        logging.info("Rename")
        path_normal_matrix_segm_tmp = pq_WithRenamedCols(path_normal_matrix_segm)
        logging.info("Load")
        normal_matrix = spark.read.parquet(path_normal_matrix_segm_tmp)
    else:
        normal_matrix = spark.read.parquet(path_normal_matrix_segm)

    # Drop junction column
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
    logging.info("partitions: {}".format(normal_matrix.rdd.getNumPartitions()))
    not_null = ' OR '.join(
        ['({} != 0.0)'.format(col_name)
         for col_name in normal_matrix.schema.names if col_name != index_name])  # SQL style  # All zeros
    normal_matrix = normal_matrix.filter(not_null)

    # Make unique
    logging.info("Make unique")
    logging.info("partitions: {}".format(normal_matrix.rdd.getNumPartitions()))
    exprs = [sf.max(sf.col(name_)).alias(name_) for name_ in normal_matrix.schema.names if name_ != index_name]
    normal_matrix = normal_matrix.groupBy(index_name).agg(*exprs)
    logging.info("partitions: {}".format(normal_matrix.rdd.getNumPartitions()))

    return spark, normal_matrix


def outlier_filtering(spark, normal_matrix, index_name, libsize_n, expr_high_limit_normal):
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
        highly_expressed_normals = ' AND '.join(['({} > {})'.format(col_name, expr_high_limit_normal)
                                                 for col_name in normal_matrix.schema.names if
                                                 col_name != index_name])  # SQL style  # Expressed kmers

        ambigous_expression_normals = ' OR '.join(['({} <= {})'.format(col_name, expr_high_limit_normal)
                                                   for col_name in normal_matrix.schema.names if
                                                   col_name != index_name])  # SQL style

    high_expr_normals = normal_matrix.filter(highly_expressed_normals).select(sf.col(index_name))
    normal_matrix = normal_matrix.filter(ambigous_expression_normals)  # TODO add condition empty matrix
    return spark, high_expr_normals, normal_matrix


def hard_filter_normals(spark, normal_matrix, index_name, libsize_n, expr_limit_normal, expr_n_limit ):
    ''' Filter normal samples based on j reads in at least n samples. The expressions are normalized for library size

    Parameters:
    ----------
    spark session
    normal matrix
    index column name
    libsize_n matrix
    expr_limit_normal (j reads)
    expr_n_limit (n samples)
    Returns :
    ----------
    spark context
    Filtered normal matrix
        '''

    if libsize_n is not None:
        normal_matrix = normal_matrix.select(index_name, *[
            sf.when((sf.col(name_) / libsize_n.loc[name_, "libsize_75percent"]) > expr_limit_normal, 1).otherwise(
                0).alias(name_)
            for name_ in normal_matrix.schema.names if name_ != index_name])
    else:
        normal_matrix = normal_matrix.select(index_name, *[
            sf.when(sf.col(name_) > expr_limit_normal, 1).otherwise(0).alias(name_)
            for name_ in normal_matrix.schema.names if name_ != index_name])

    logging.info("partitions: {}".format(normal_matrix.rdd.getNumPartitions()))
    logging.info("reduce")
    rowsum = (reduce(add, (sf.col(x) for x in normal_matrix.columns[1:]))).alias("sum")
    logging.info("partitions: {}".format(normal_matrix.rdd.getNumPartitions()))
    logging.info("filter")
    normal_matrix = normal_matrix.filter(rowsum >= expr_n_limit)

    return spark, normal_matrix


def process_cancers(spark, cancer_path, cancer_sample, drop_cols, expression_fields, jct_col, index_name, libsize_c, cross_junction):
    ''' Preprocess cancer samples
    - Make kmers unique
    - Filter kmers on junction status
    - Normalize
    Parameters:
    ----------
    spark: spark session
    cancer path: path to cancer file
    cancer_sample: associated cacncer ID
    drop_cols: colums to be dropped
    expression_fields: list of segment and junction expression column names
    jct_col: junction status column name
    index_name: kmer column name
    libsize_c: libsize matrix for cancer samples
    cross_junction: Information to filter on juction status. None (both, no filtering), True (junction), False (non junction)
    Returns
    --------
    spark: spark session
    cancer_kmers: cancer kmers matrix,
    cancer_path_tmp: path of renamed temporary file
    jct_type: string indicating which junction filtering has been performed
    '''

    def collapse_values(value):
        return max([np.float(i) if i != 'nan' else 0.0 for i in value.split('/')])  # np.nanmax not supported

    cancer_path_tmp = pq_WithRenamedCols(cancer_path)
    cancer_kmers = spark.read.parquet(cancer_path_tmp)
    logging.info("Renamed done")

    # Drop junction column
    for drop_col in drop_cols:
        cancer_kmers = cancer_kmers.drop(sf.col(drop_col))
    logging.info("Collapse kmer horizontal")

    # Remove the '/' in the expression data (kmers duplicate within a gene have 'expression1/expression2' format
    local_max = sf.udf(collapse_values, st.FloatType())
    for name_ in expression_fields:
        cancer_kmers = cancer_kmers.withColumn(name_, local_max(name_))

    # Filter on juction status
    if cross_junction == 1:
        logging.info("Keep cross junction kmers")
        cancer_kmers = cancer_kmers.filter("{} == True".format(jct_col))
        jct_type = "only-jct"
    elif cross_junction == 0:
        logging.info("Keep non-cross junction kmers")
        cancer_kmers = cancer_kmers.filter("{} == False".format(jct_col))
        jct_type = "non-jct"
    else:
        logging.info("Keep cross and non- junction kmers")
        jct_type = "no-jct-filt"
    if len(cancer_kmers.head(1)) == 0:
        logging.info("WARNING: Restriction on junction status removed all foreground... Exiting")
        sys.exit(1)

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

    return spark, cancer_kmers, cancer_path_tmp, jct_type


def remove_uniprot(spark, cancer_kmers, uniprot, index_name):
    def I_L_replace(value):
        return value.replace('I', 'L')
    if uniprot is not None:
        uniprot = spark.read.csv(uniprot, sep='\t', header=None)
        uniprot_header = index_name + "_IL_eq"
        uniprot = uniprot.withColumnRenamed("_c0", uniprot_header)
        # Make isoleucine and leucine equivalent
        I_L_equ = sf.udf(I_L_replace, st.StringType())
        uniprot = uniprot.withColumn(uniprot_header, I_L_equ(uniprot_header))
        cancer_kmers = cancer_kmers.withColumn(uniprot_header, I_L_equ(index_name))
        cancer_kmers = cancer_kmers.join(uniprot, cancer_kmers[uniprot_header] == uniprot[uniprot_header],
                                         how='left_anti')
    return  spark, cancer_kmers


def save_spark(cancer_kmers, output_dir, path_final_fil):
    # save
    if len(cancer_kmers.head(1)) > 0:
        logging.info("Save to {}".format(path_final_fil))
        pathlib.Path(output_dir).mkdir(exist_ok=True, parents=True)

        cancer_kmers.repartition(1).write.mode('overwrite').parquet(path_final_fil)
    else:
        logging.info("WARNING: no saving performed, {} would be empty".format(path_final_fil))





### Main
def mode_cancerspecif(arg):

    spark_cfg = default_spark_config(arg.cores, arg.mem_per_core)
    with create_spark_session_from_config(spark_cfg) as spark:
        #if os.path.exists(os.path.join(arg.output_dir, "checkpoint")):
        #    shutil.rmtree(os.path.join(arg.output_dir, "checkpoint"))
        #pathlib.Path(os.path.join(arg.output_dir, "checkpoint")).mkdir(exist_ok=True, parents=True)
        #spark.sparkContext.setCheckpointDir(os.path.join(arg.output_dir, "checkpoint"))

        ### Preprocessing Normals
        logging.info("\n >>>>>>>> Preprocessing Normal samples")
        index_name = 'kmer'
        jct_col = "iscrossjunction"
        spark, normal_matrix = process_normals(spark, index_name, jct_col, arg.path_normal_matrix_segm, arg.whitelist)

        ### Preprocessing Libsize
        logging.info("\n >>>>>>>> Preprocessing libsizes")
        if arg.path_normal_libsize:
            libsize_n = process_libsize(arg.path_normal_libsize)
        else:
            libsize_n = None
        if arg.path_cancer_libsize:
            libsize_c = process_libsize(arg.path_cancer_libsize)


        ### NORMALS: Statistical Filtering
        if arg.statistical:
            logging.info("\n >>>>>>>> Normals: Perform statistical filtering")
            # Remove outlier kmers before statistical modelling (Very highly expressed kmers do not follow a NB, we classify them as expressed without hypothesis testing)
            if arg.expr_high_limit_normal is not None:
                logging.info( "... Normal kmers with expression >= {} in >= 1 sample are truly expressed. \n Will be substracted from cancer set".format(arg.expr_high_limit_normal))
                spark, high_expr_normals, normal_matrix = outlier_filtering(spark, normal_matrix, index_name, libsize_n, arg.expr_high_limit_normal)

            if arg.tissue_grp_files is not None:
                modelling_grps = []
                for tissue_grp in arg.tissue_grp_files:
                    grp = pd.read_csv(tissue_grp, header = None)[0].to_list()
                    grp = [name_.replace('-','').replace('.','').replace('_','') for name_ in grp]
                    grp.append(index_name)
                    modelling_grps.append(grp)
            else:
                modelling_grps = [[name_ for name_ in normal_matrix.schema.names if name_ != index_name]]

            logging.info(">>>... Fit Negative Binomial distribution on normal kmers ")
            for grp in modelling_grps:
            # Fit NB and Perform hypothesis testing
                spark, normal_matrix = fit_NB(spark, normal_matrix, index_name, arg.output_dir, arg.path_normal_matrix_segm, libsize_n, arg.cores)
                normal_matrix = normal_matrix.filter(sf.col("proba_zero") < arg.threshold_noise) # Expressed kmers

            # Join on the kmers segments. Take the kmer which junction expression is not zero everywhere

        ### NORMALS: Hard Filtering
        else:
            logging.info("\n >>>>>>>> Normals: Perform Hard Filtering \n (expressed in {} samples with {} normalized counts)".format(arg.expr_n_limit, arg.expr_limit_normal))
            logging.info("expression filter")
            spark, normal_matrix = hard_filter_normals(spark, normal_matrix, index_name, libsize_n, arg.expr_limit_normal, arg.expr_n_limit)

            save_filtered_normals = True
            if save_filtered_normals:
                extension = '.pq'
                path_tmp_n = os.path.join(arg.output_dir, os.path.basename(arg.path_normal_matrix_segm).split('.')[
                    0] + '_expr-in-{}-samples-with-{}-normalized-cts'.format(arg.expr_n_limit, arg.expr_limit_normal) + extension)
                save_spark(normal_matrix, arg.output_dir, path_tmp_n)
                logging.info(path_tmp_n)
            normal_matrix = normal_matrix.select(sf.col(index_name))

        ### Apply filtering to foreground
        if arg.expression_fields_c is None:
            expression_fields_orig =  ['segment_expr', 'junction_expr']
        else:
            expression_fields_orig = arg.expression_fields_c

        expression_fields = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in expression_fields_orig]
        drop_cols = ['id']

        for cancer_path, cancer_sample in zip(arg.paths_cancer_samples, arg.ids_cancer_samples):
            cancer_sample = cancer_sample.replace('-', '').replace('.', '').replace('_', '')

            # Preprocess cancer samples
            logging.info("\n >>>>>>>> Cancers: Perform differential filtering sample {}".format(cancer_sample))
            spark, cancer_kmers, cancer_path_tmp, jct_type = process_cancers(spark, cancer_path, cancer_sample, drop_cols,
                                                                      expression_fields, jct_col, index_name,
                                                                        libsize_c, arg.cross_junction)

            ###Perform Background removal
            #logging.info("Perform join")
            #logging.info("BROADCAST NORMALS")
            #normal_matrix = sf.broadcast(normal_matrix)
            #cancer_kmers = sf.broadcast(cancer_kmers)
            #cancer_kmers = cancer_kmers.join(normal_matrix, cancer_kmers["kmer"] == normal_matrix["kmer"], how='left_anti')

            extension = '.pq'
            path_tmp_c = os.path.join(arg.output_dir, os.path.basename(arg.paths_cancer_samples[0]).split('.')[
                0] + '_expressed_normalized_' + jct_type + extension)
            save_spark(cancer_kmers, arg.output_dir, path_tmp_c)
            logging.info(path_tmp_c)

            logging.info("Filtering normal background")
            cancer_kmers = cancer_kmers.join(normal_matrix, cancer_kmers["kmer"] == normal_matrix["kmer"], how='left_anti')

            path_final_fil = os.path.join(arg.output_dir, os.path.basename(arg.paths_cancer_samples[0]).split('.')[
                0] + '_'+ jct_type + '_filter-for-normals'  + extension)
            save_spark(cancer_kmers, arg.output_dir, path_final_fil)


            ### Remove Uniprot
            logging.info("Filtering kmers in uniprot")
            spark, cancer_kmers = remove_uniprot(spark, cancer_kmers, arg.uniprot, index_name)
            path_final_fil = os.path.join(arg.output_dir, os.path.basename(arg.paths_cancer_samples[0]).split('.')[
                0] + '_'+ jct_type + '_filter-for-normals-and-uniprot' + extension)
            save_spark(cancer_kmers, arg.output_dir, path_final_fil)


            os.remove(cancer_path_tmp)




            #TODO Implement the intersection of the modelling tissues
            #TODO transfer to junction xpression

