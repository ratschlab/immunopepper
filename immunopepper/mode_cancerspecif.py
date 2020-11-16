import logging
import numpy as np
import os
import pandas as pd
from pyspark.sql import functions as sf
from pyspark.sql import types as st
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, Formula, r
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
import scipy
from scipy import stats


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

    #count_matrix.rownames = count_matrix.get_attrib('kmer')
    design_formula = Formula(' ~ 1')


    dds0 = deseq.DESeqDataSetFromMatrix(countData=count_matrix,
                                            colData=design_matrix,
                                            design=design_formula)
    dds0 = BiocGenerics.estimateSizeFactors(dds0, type="poscounts")
    order_size_factor = list(dds0.do_slot('colData').do_slot('rownames'))
    dds0.do_slot('colData').do_slot('listData')[1] = ro.vectors.FloatVector(list(normalize.loc[order_size_factor, 'libsize_75percent'])) # Enforce size factors



    dds = deseq.DESeq(dds0, parallel=True, BPPARAM=BiocParallel.MulticoreParam(cores),
                      sfType="poscounts", # Will run 1. estimation of size factors: estimateSizeFactors # parameter "poscounts"
                      fitType="parametric" # 2. estimation of dispersion: estimateDispersions # parameter "parametric"
                      )

    #normalized_count_matrix = deseq.counts(self, normalized=Tru os.path.join(cancer_dir, 'expression_counts.libsize.tsv'),e)
    deseq_result = deseq.results(dds)
    fit_res = to_dataframe(deseq_result)
    disp = to_dataframe(deseq.dispersions(dds)).rename({'x': 'dispersion'}, axis = 1)
    disp.index = fit_res.index
    fit_res = pd.concat([fit_res['baseMean'], disp], axis=1)


    return fit_res



def filter_cancer_samples(cancer_sample_paths, normal_kmers_expressed, spark_session):
    for cancer_sample in cancer_sample_paths:
        normal_kmers = spark_session.read.parquet(cancer_sample)
        cancer_sample.join(normal_kmers_expressed.select('kmer'), how='inner', on=['kmer'])  # Remark left_anti could be useful



def mode_cancerspecif(arg):
    spark_cfg = default_spark_config(arg.cores, arg.mem_per_core)

    with create_spark_session_from_config(spark_cfg) as spark:
        def nb_cdf(mean_, disp):
            probA = disp / (disp + mean_)
            N = (probA * mean_) / (1 - probA)
            return float(scipy.stats.nbinom.cdf(0, n=N, p=probA))

        #spark.sparkContext.addFile('/Users/laurieprelot/software/anaconda/lib/python3.6/site-packages/scipy.zip') # TODO do we need it ? packages in config
        ### With Pandas
        # normal_matrix = pd.read_parquet(arg.path_normal_matrix_segm, engine = 'pyarrow') #
        # normal_matrix  = normal_matrix.drop(['is_cross_junction'], axis=1)
        # normal_matrix = normal_matrix.set_index('kmer')
        # normal_matrix = normal_matrix.astype(float)
        # normal_matrix[normal_matrix.isnull()] = 0.0
        # normal_matrix.columns = ['TCGA-AO-A12D-01A-11', 'TCGA-AR-A0TT-01A-31'] #TODO remove
        # normal_matrix = normal_matrix.reset_index('kmer')


        ### Preprocessing Normals
        index_name = 'kmer'
        normal_matrix = spark.read.parquet(arg.path_normal_matrix_segm)
        normal_matrix = normal_matrix.drop('is_cross_junction')
        for name_ in normal_matrix.schema.names:
            if name_ != index_name:
                normal_matrix = normal_matrix.withColumn(name_, sf.col(name_).cast(st.DoubleType()))
                normal_matrix = normal_matrix.na.fill(0)

        ### Preprocessing Libsize
        libsize = pd.read_csv(arg.path_normal_libsize, sep = '\t')
        libsize['libsize_75percent'] = libsize['libsize_75percent'] / np.median(libsize['libsize_75percent'])
        libsize['sample'] = [sample.replace('-', '_') for sample in libsize['sample']]
        libsize = libsize.set_index('sample')



        if arg.statistical: ### Statistical Filtering

            ### Remove outlier kmers before statistical modelling (Very highly expressed kmers do not follow a NB, we classify them as expressed without hypothesis testing)
            if arg.expr_high_limit_normal is not None:
                logging.info( ">>>>>>>> Normal kmers with expression >= {} in >= 1 sample are truly expressed. Will be substracted from cancer set".format(arg.expr_high_limit_normal))
                # normal_matrix = spark.createDataFrame(normal_matrix)
                for name_ in normal_matrix.schema.names:
                    normal_matrix = normal_matrix.withColumnRenamed(name_, name_.replace('-',
                                                                                             '_'))  # '-' causing issue in filter function
                highly_expressed_normals = ' AND '.join(
                    ['({} > {})'.format(col_name, arg.expr_high_limit_normal) for col_name in
                     normal_matrix.schema.names if col_name != index_name])  # SQL style  # Expressed kmers
                ambigous_expression_normals = ' OR '.join(
                    ['({} <= {})'.format(col_name, arg.expr_high_limit_normal) for col_name in
                     normal_matrix.schema.names if col_name != index_name])  # SQL style
                high_expr_normals = normal_matrix.filter(highly_expressed_normals).select(sf.col(index_name))
                normal_matrix = normal_matrix.filter(ambigous_expression_normals)  # TODO add condition empty matrix

            ### Perform hypothesis testing
            logging.info(">>>>>>>> Fit Negative Binomial distribution on normal kmers ")
            design_matrix = pd.DataFrame([1] * (len(normal_matrix.columns) - 1), columns=["design"])
            design_matrix['sample'] = [col_name for col_name in normal_matrix.schema.names if col_name != index_name]
            design_matrix = design_matrix.set_index('sample')
            res = DESeq2(normal_matrix.toPandas().set_index(index_name), design_matrix, normalize=libsize, cores=1)
            save_pd_toparquet(os.path.join(arg.output_dir, os.path.basename(arg.path_normal_matrix_segm).split('.')[0] + 'deseq_fit' + '.pq.gz'), res, compression = 'gzip', verbose = False)
            # res = pd.read_parquet(os.path.join(arg.output_dir, os.path.basename(arg.path_path_normal_matrix_segm).split('.')[0] + 'deseq_fit' + '.pq.gz'))
            res_2 = spark.createDataFrame(res.reset_index())
            stat_udf = sf.udf(nb_cdf, st.FloatType())
            res_2 = res_2.withColumn("proba_zero", stat_udf(sf.col('baseMean'), sf.col('dispersion')))
            expressed_normals_dsq = res_2.filter(sf.col("proba_zero") < arg.threshold_noise) # Expressed kmers
            #TODO transfer to junction xpression

            ### Apply filtering to foreground
            filter_cancer_samples(arg.paths_cancer_samples , expressed_normals_dsq, spark)


        else: ### Hard Filtering
            for name_ in normal_matrix.schema.names:
                if name_ != index_name:
                    normal_matrix = normal_matrix.withColumn(name_, sf.when(sf.col(name_) > arg.expr_limit_normal, 1).otherwise(0))
            normal_matrix  = normal_matrix.withColumn('passing_expr_criteria', sum(normal_matrix[name_] for name_ in normal_matrix.schema.names if  name_ != index_name))
            normal_matrix = normal_matrix.filter(sf.col("passing_expr_criteria") >= arg.expr_n_limit)
            ### Apply filtering to foreground 
            filter_cancer_samples(arg.paths_cancer_samples, normal_matrix, spark)




