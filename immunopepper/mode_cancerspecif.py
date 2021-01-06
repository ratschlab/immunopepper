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
    #count_matrix.rownames = count_matrix.get_attrib('kmer')
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
    #normalized_count_matrix = deseq.counts(self, normalized=Tru os.path.join(cancer_dir, 'expression_counts.libsize.tsv'),e)

    deseq_result = deseq.results(dds)
    fit_res = to_dataframe(deseq_result)
    disp = to_dataframe(deseq.dispersions(dds)).rename({'x': 'dispersion'}, axis = 1)
    disp.index = fit_res.index
    fit_res = pd.concat([fit_res['baseMean'], disp], axis=1)
    return fit_res


def preprocess_libsize(path_lib):
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
    logging.info("resaving to {}".format(path_tmp))
    pqwriter = pq.ParquetWriter(path_tmp, df.schema, compression=None)
    pqwriter.write_table(df)
    pqwriter.close()
    return path_tmp

def mode_cancerspecif(arg):
    spark_cfg = default_spark_config(arg.cores, arg.mem_per_core)

    with create_spark_session_from_config(spark_cfg) as spark:
        if os.path.exists(os.path.join(arg.output_dir, "checkpoint")):
            shutil.rmtree(os.path.join(arg.output_dir, "checkpoint"))
        pathlib.Path(os.path.join(arg.output_dir, "checkpoint")).mkdir(exist_ok=True, parents=True)
        spark.sparkContext.setCheckpointDir(os.path.join(arg.output_dir, "checkpoint"))

        def nb_cdf(mean_, disp):
            probA = disp / (disp + mean_)
            N = (probA * mean_) / (1 - probA)
            return float(scipy.stats.nbinom.cdf(0, n=N, p=probA))

        def collapse_values(value):
            return  max([np.float(i) if i!='nan' else 0.0 for i in value.split('/') ]) #np.nanmax not supported

        def cast_type_dbl(normal_matrix, name_list, index_name ):
            return normal_matrix.select( [sf.col(name_).cast(st.DoubleType()).alias(name_) if name_ != index_name else sf.col(name_) for name_ in name_list] )

        #spark.sparkContext.addFile('/Users/laurieprelot/software/anaconda/lib/python3.6/sitesddpackages/scipy.zip') # TODO do we need it ? packages in config

        ### Preprocessing Normals
        logging.info(">>>>>>>> Preprocessing Normal samples")
        index_name = 'kmer'
        jct_col = "iscrossjunction"
        # Rename
        rename = False #For development 
        if rename:
            logging.info("Rename")
            path_normal_matrix_segm_tmp =  pq_WithRenamedCols(arg.path_normal_matrix_segm)
            logging.info("Load")
            normal_matrix = spark.read.parquet(path_normal_matrix_segm_tmp)
        else: 
            normal_matrix = spark.read.parquet(arg.path_normal_matrix_segm)


        #logging.info("repartition")
        #normal_matrix =  normal_matrix.repartition(arg.cores * 10)
        #logging.info("partitions: {}".format(normal_matrix.rdd.getNumPartitions()))

        normal_matrix = normal_matrix.drop(jct_col)
        
        #logging.info("cache")
        #normal_matrix.cache()
        #logging.info("checkpoint")
        #normal_matrix.checkpoint()

        # Cast type and fill nans
        logging.info("Cast types")

        if arg.whitelist is not None:
            whitelist = pd.read_csv(arg.whitelist, sep = '\t', header = None)[0].to_list()
            whitelist = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in whitelist]
            whitelist.append(index_name)
            normal_matrix = cast_type_dbl(normal_matrix, whitelist, index_name)
        else:
            normal_matrix =  cast_type_dbl(normal_matrix, normal_matrix.schema.names, index_name)


        logging.info("Remove Nans")
        normal_matrix = normal_matrix.na.fill(0)

        logging.info("Remove non expressed kmers SQL-")

        #logging.info("cache")
        #normal_matrix.cache()
        #logging.info("checkpoint")
        #normal_matrix.checkpoint()
        logging.info("partitions: {}".format(normal_matrix.rdd.getNumPartitions()))


        # Remove lines of Nans (Only present in junction file)
        # normal_matrix = normal_matrix.withColumn('allnull', sum(normal_matrix[name_] for name_ in normal_matrix.schema.names if name_ != index_name))
        # normal_matrix = normal_matrix.filter(sf.col("allnull") > 0.0 )
        # normal_matrix = normal_matrix.drop("allnull")

        not_null = ' OR '.join(
            ['({} != 0.0)'.format(col_name)
        for col_name in normal_matrix.schema.names if col_name != index_name])  # SQL style  # All zeros
        normal_matrix = normal_matrix.filter(not_null)

        #logging.info("cache")
        #normal_matrix.cache()
        #logging.info("checkpoint")
        #normal_matrix.checkpoint()
        logging.info("partitions: {}".format(normal_matrix.rdd.getNumPartitions()))

        logging.info("Make unique")
        # Make unique

        exprs = [sf.max(sf.col(name_)).alias(name_) for name_ in  normal_matrix.schema.names if name_ != index_name]
        normal_matrix = normal_matrix.groupBy(index_name).agg(*exprs)

        #logging.info("cache")
        #normal_matrix.cache()
        #logging.info("checkpoint")
        #normal_matrix.checkpoint()
        logging.info("partitions: {}".format(normal_matrix.rdd.getNumPartitions()))

        ### Preprocessing Libsize
        logging.info(">>>>>>>> Preprocessing libsizes")
        if arg.path_normal_libsize:
            libsize_n = preprocess_libsize(arg.path_normal_libsize)
        else:
            libsize_n = None
        if arg.path_cancer_libsize:
            libsize_c = preprocess_libsize(arg.path_cancer_libsize)


        ### NORMALS: Statistical Filtering
        if arg.statistical:
            logging.info(">>>>>>>>  Perform statistical filtering")
            # Remove outlier kmers before statistical modelling (Very highly expressed kmers do not follow a NB, we classify them as expressed without hypothesis testing)
            if arg.expr_high_limit_normal is not None:
                logging.info( "... Normal kmers with expression >= {} in >= 1 sample are truly expressed. \n Will be substracted from cancer set".format(arg.expr_high_limit_normal))
                # normal_matrix = spark.createDataFrame(normal_matrix)
                if libsize_n is not None:
                    highly_expressed_normals = ' AND '.join( ['({} > {})'.format(col_name, arg.expr_high_limit_normal * libsize_n.loc[col_name, "libsize_75percent"])
                         for col_name in normal_matrix.schema.names if col_name != index_name])  # SQL style  # Expressed kmers

                    ambigous_expression_normals = ' OR '.join(['({} <= {})'.format(col_name, arg.expr_high_limit_normal  * libsize_n.loc[col_name, "libsize_75percent"])
                         for col_name in normal_matrix.schema.names if col_name != index_name])  # SQL style
                else:
                    highly_expressed_normals = ' AND '.join( ['({} > {})'.format(col_name, arg.expr_high_limit_normal )
                         for col_name in normal_matrix.schema.names if col_name != index_name])  # SQL style  # Expressed kmers

                    ambigous_expression_normals = ' OR '.join(['({} <= {})'.format(col_name, arg.expr_high_limit_normal )
                         for col_name in normal_matrix.schema.names if col_name != index_name])  # SQL style

                high_expr_normals = normal_matrix.filter(highly_expressed_normals).select(sf.col(index_name))
                normal_matrix = normal_matrix.filter(ambigous_expression_normals)  # TODO add condition empty matrix

                #normal_matrix.cache()
                #normal_matrix.checkpoint()

            if arg.tissue_grp_files is not None:
                modelling_grps = []
                for tissue_grp in arg.tissue_grp_files:
                    grp = pd.read_csv(tissue_grp, header = None)[0].to_list()
                    grp = [name_.replace('-','').replace('.','').replace('_','') for name_ in grp]
                    grp.append(index_name)
                    modelling_grps.append(grp)

            else:
                modelling_grps = [[name_ for name_ in normal_matrix.schema.names if name_ != index_name]]

            for grp in modelling_grps:
            # Perform hypothesis testing
                logging.info(">>>... Fit Negative Binomial distribution on normal kmers ")
                design_matrix = pd.DataFrame([1] * (len(normal_matrix.columns) - 1), columns=["design"])
                design_matrix['sample'] = [col_name for col_name in normal_matrix.schema.names if col_name != index_name]

                # Run DESEq2
                design_matrix = design_matrix.set_index('sample')
                logging.info("Run DESeq2")
                normal_matrix = DESeq2(normal_matrix.toPandas().set_index(index_name), design_matrix, normalize=libsize_n, cores=arg.cores)
                save_pd_toparquet(os.path.join(arg.output_dir, os.path.basename(arg.path_normal_matrix_segm).split('.')[0] + 'deseq_fit' + '.pq.gz'), normal_matrix, compression = 'gzip', verbose = False)
                # Test on probability of noise
                logging.info("Test if noise")
                normal_matrix = spark.createDataFrame(normal_matrix.reset_index())
                stat_udf = sf.udf(nb_cdf, st.FloatType())
                logging.info("Filter out noise from normals")
                normal_matrix = normal_matrix.withColumn("proba_zero", stat_udf(sf.col('baseMean'), sf.col('dispersion')))
                normal_matrix = normal_matrix.filter(sf.col("proba_zero") < arg.threshold_noise) # Expressed kmers

            # Join on the kmers segments. Take the kmer which junction expression is not zero everywhere

        ### NORMALS: Hard Filtering
        else:
            logging.info(">>>>>>>> Hard Filtering")
        # Filter based on expression level in at least n samples

            logging.info("expression filter")

            if libsize_n is not None:
                normal_matrix= normal_matrix.select(index_name, *[
                    sf.when((sf.col(name_) / libsize_n.loc[name_, "libsize_75percent"]) > arg.expr_limit_normal, 1).otherwise(0).alias(name_)
                    for name_ in normal_matrix.schema.names if name_ != index_name])
            else:
                normal_matrix= normal_matrix.select(index_name, *[
                    sf.when(sf.col(name_)  > arg.expr_limit_normal, 1).otherwise(0).alias(name_)
                    for name_ in normal_matrix.schema.names if name_ != index_name])
            
            #logging.info("cache")
            #normal_matrix.cache()
            #logging.info("checkpoint")
            #normal_matrix.checkpoint()
            logging.info("partitions: {}".format(normal_matrix.rdd.getNumPartitions()))    
            
            logging.info("reduce")
            rowsum = (reduce(add, (sf.col(x) for x in normal_matrix.columns[1:]))).alias("sum")

            #logging.info("cache")
            #normal_matrix.cache()
            #logging.info("checkpoint")
            #normal_matrix.checkpoint()
            logging.info("partitions: {}".format(normal_matrix.rdd.getNumPartitions()))
            
            logging.info("filter")
            normal_matrix = normal_matrix.filter(rowsum >= arg.expr_n_limit)

            # for name_ in normal_matrix.schema.names:
            #     if name_ != index_name:
            #         normal_matrix = normal_matrix.withColumn(name_, sf.when(
            #             sf.col(name_) / libsize_n.loc[name_, "libsize_75percent"] > arg.expr_limit_normal, 1).otherwise(
            #             0))  # Normalize on the fly
            # normal_matrix = normal_matrix.withColumn('passing_expr_criteria', sum(
            #     normal_matrix[name_] for name_ in normal_matrix.schema.names if name_ != index_name))
            # normal_matrix = normal_matrix.filter(sf.col("passing_expr_criteria") >= arg.expr_n_limit)

        #normal_matrix.cache()
        #normal_matrix.checkpoint()

        logging.info("No caching")
        ### Apply filtering to foreground
        if arg.expression_fields_c is None:
            expression_fields_orig =  ['segment_expr', 'junction_expr']
        else:
            expression_fields_orig = arg.expression_fields_c
        expression_fields = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in expression_fields_orig]
        drop_cols = ['id']
        for cancer_path, cancer_sample in zip(arg.paths_cancer_samples, arg.ids_cancer_samples):
            cancer_sample = cancer_sample.replace('-', '').replace('.', '').replace('_', '')
            logging.info("... Process cancer sample {}".format(cancer_sample))
            cancer_path_tmp = pq_WithRenamedCols(cancer_path)
            cancer_kmers = spark.read.parquet(cancer_path_tmp)
            logging.info("Renamed done")
            for drop_col in drop_cols:
                cancer_kmers = cancer_kmers.drop(sf.col(drop_col))
            logging.info("Collapse kmer horizontal")
            # Remove the '/' in the data
            local_max = sf.udf(collapse_values, st.FloatType())
            for name_ in expression_fields:
                cancer_kmers = cancer_kmers.withColumn(name_, local_max(name_))

            if arg.cross_junction == 1:
                logging.info("Keep cross junction kmers")
                cancer_kmers = cancer_kmers.filter("{} == True".format(jct_col))
            elif arg.cross_junction == 0:
                logging.info("Keep non-cross junction kmers")
                cancer_kmers = cancer_kmers.filter("{} == False".format(jct_col))
            else:
                logging.info("Keep cross and non- junction kmers")

            if len(cancer_kmers.head(1)) == 0:
                logging.info("WARNING: Restriction on junction status removed all foreground... Exiting")
                sys.exit(1)


            logging.info("Collapse kmer vertical")
            # Make unique, max
            cancer_kmers = cancer_kmers.withColumn(jct_col, sf.col(jct_col).cast("boolean").cast("int"))
            exprs = [sf.max(sf.col(name_)).alias(name_) for name_ in cancer_kmers.schema.names if name_ != index_name]
            cancer_kmers = cancer_kmers.groupBy(index_name).agg(*exprs)
            #cancer_kmers = cancer_kmers.sort(sf.col(index_name))

            #Remove double zeros filelds
            cancer_kmers = cancer_kmers.withColumn('allnull', sum(cancer_kmers[name_] for name_ in expression_fields))
            cancer_kmers = cancer_kmers.filter(sf.col("allnull") > 0.0)
            cancer_kmers = cancer_kmers.drop("allnull")

            #Normalize by library size
            if libsize_c is not None:
                cancer_kmers = cancer_kmers.withColumn(name_, sf.round(cancer_kmers[name_] /libsize_c.loc[cancer_sample, "libsize_75percent"], 2))
            else:
                cancer_kmers = cancer_kmers.withColumn(name_, sf.round( cancer_kmers[name_] , 2))

            #Perform Background removal
            logging.info("Perform join")
            logging.info("NO BROADCAST")
            #normal_matrix = sf.broadcast(normal_matrix)
            #cancer_kmers = sf.broadcast(cancer_kmers)
            cancer_kmers = cancer_kmers.join(normal_matrix, cancer_kmers["kmer"] == normal_matrix["kmer"], how='left_anti')

            if len(cancer_kmers.head(0)) > 0:
                logging.info("save filtered output")
                pathlib.Path(arg.output_dir).mkdir(exist_ok= True, parents= True)
                extension = '.pq'
                path_final_fil = os.path.join(arg.output_dir, os.path.basename(arg.paths_cancer_samples[0]).split('.')[0] + '_no_normals' + extension)
                cancer_kmers.repartition(1).write.mode('overwrite').parquet(path_final_fil)
            else:
                logging.info("WARNING: filtering removed all Foreground")
            os.remove(cancer_path_tmp)


            #TODO Implement the intersection of the modelling tissues
            #TODO Implement intermediary saving # no collect slow
            #TODO transfer to junction xpression
            #TODO Uniprot filtering

