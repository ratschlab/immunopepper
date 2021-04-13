import logging
import os
import pandas as pd
from pyspark.sql import functions as sf


from .config import create_spark_session_from_config
from .config import default_spark_config
from .statistical import process_normals
from .statistical import save_spark
from .statistical import outlier_filtering
from .statistical import hard_filter_normals
from .statistical import fit_NB
from .statistical import pq_WithRenamedCols
from .statistical import filter_cancer
from .statistical import preprocess_cancers
from .statistical import remove_uniprot
from .statistical import process_libsize
from .statistical import filter_statistical


### Main
def mode_cancerspecif(arg):

    spark_cfg = default_spark_config(arg.cores, arg.mem_per_core, arg.parallelism, tmp_dir=arg.scratch_dir)
    with create_spark_session_from_config(spark_cfg) as spark:
        # if os.path.exists(os.path.join(arg.output_dir, "checkpoint")):
        #     shutil.rmtree(os.path.join(arg.output_dir, "checkpoint"))
        # pathlib.Path(os.path.join(arg.output_dir, "checkpoint")).mkdir(exist_ok=True, parents=True)
        # spark.sparkContext.setCheckpointDir(os.path.join(arg.output_dir, "checkpoint"))

        index_name = 'kmer'
        jct_col = "iscrossjunction"
        save_intermed = False
        save_kmersnormal = False
        save_canc_int = False

        ### Preprocessing Libsize
        logging.info("\n >>>>>>>> Preprocessing libsizes")
        if arg.path_normal_libsize:
            libsize_n = process_libsize(arg.path_normal_libsize)
        else:
            libsize_n = None
        if arg.path_cancer_libsize:
            libsize_c = process_libsize(arg.path_cancer_libsize)

        if arg.path_normal_kmer_list is None:
            ### Preprocessing Normals
            logging.info("\n >>>>>>>> Preprocessing Normal samples")

            normal_matrix = process_normals(spark, index_name, jct_col, arg.path_normal_matrix_segm, arg.output_dir, arg.whitelist, arg.parallelism, cross_junction = 0).union(process_normals(spark, index_name, jct_col, arg.path_normal_matrix_edge, arg.output_dir, arg.whitelist, arg.parallelism, cross_junction = 1))
            if save_intermed:
                path_ = os.path.join(arg.output_dir, 'normals_merge-segm-edge.tsv')
                save_spark(normal_matrix, arg.output_dir, path_)

            # Take max expression between edge or segment expression
            exprs = [sf.max(sf.col(name_)).alias(name_) for name_ in normal_matrix.schema.names if name_ != index_name]
            normal_matrix = normal_matrix.groupBy(index_name).agg(*exprs)
            if save_intermed:
                path_ = os.path.join(arg.output_dir, 'normals_merge-segm-edge_max_uniq.tsv')
                save_spark(normal_matrix, arg.output_dir, path_)


            ### NORMALS: Statistical Filtering
            if arg.statistical:
                logging.info("\n >>>>>>>> Normals: Perform statistical filtering")
                # Remove outlier kmers before statistical modelling (Very highly expressed kmers do not follow a NB, we classify them as expressed without hypothesis testing)
                if arg.expr_high_limit_normal is not None:
                    logging.info( "... Normal kmers with expression >= {} in >= 1 sample are truly expressed. \n Will be substracted from cancer set".format(arg.expr_high_limit_normal))
                    high_expr_normals, normal_matrix = outlier_filtering(spark, normal_matrix, index_name, libsize_n, arg.expr_high_limit_normal)

                filter_statistical(spark, arg.tissue_grp_files, normal_matrix, index_name, arg.path_normal_matrix_segm,
                                       libsize_n, arg.threshold_noise, arg.output_dir, arg.cores)
            ### NORMALS: Hard Filtering
            else:
                logging.info("\n >>>>>>>> Normals: Perform Hard Filtering \n (expressed in {} samples with {} normalized counts)".format(arg.n_samples_lim_normal, arg.expr_limit_normal))
                logging.info("expression filter")

                path_normal_kmers = hard_filter_normals(normal_matrix, index_name, libsize_n, arg.output_dir, arg.expr_limit_normal, arg.n_samples_lim_normal)
                normal_matrix = spark.read.csv(path_normal_kmers, sep=r'\t', header=False)
        else:
            normal_matrix = spark.read.csv(arg.path_normal_kmer_list, sep=r'\t', header=False)
        normal_matrix = normal_matrix.withColumnRenamed('_c0', index_name)
        normal_matrix = normal_matrix.select(sf.col(index_name))


        ### Apply filtering to foreground
        if arg.expression_fields_c is None:
            expression_fields_orig =  ['segment_expr', 'junction_expr'] # segment 1 , junction 2
        else:
            expression_fields_orig = arg.expression_fields_c

        expression_fields = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in expression_fields_orig]
        drop_cols = ['id']

        for cancer_path, cancer_sample in zip(arg.paths_cancer_samples, arg.ids_cancer_samples):
            cancer_sample = cancer_sample.replace('-', '').replace('.', '').replace('_', '')
            cancer_path_tmp = pq_WithRenamedCols(cancer_path, arg.output_dir)
            cancer_kmers = spark.read.parquet(cancer_path_tmp)

            logging.info("Renamed done")

            # Preprocess cancer samples
            logging.info("\n >>>>>>>> Cancers: Perform differential filtering sample {}".format(cancer_sample))

            cancer_kmers = filter_cancer(preprocess_cancers(cancer_kmers, cancer_sample, drop_cols,
                                                                      expression_fields, jct_col, index_name,
                                                                        libsize_c, 1),
                                         preprocess_cancers(cancer_kmers, cancer_sample, drop_cols,
                                                            expression_fields, jct_col, index_name,
                                                            libsize_c, 0),
                                         index_name, expression_fields, arg.expr_limit_cancer, arg.parallelism)


            if save_canc_int:
                extension = '.tsv'
                path_tmp_c = os.path.join(arg.output_dir, os.path.basename(arg.paths_cancer_samples[0]).split('.')[
                0] + '_expressed_normalized_'  + extension)
                save_spark(cancer_kmers, arg.output_dir, path_tmp_c)
            
            logging.info("partitions: {}".format(cancer_kmers.rdd.getNumPartitions()))

            logging.info("Filtering normal background")
            cancer_kmers = cancer_kmers.join(normal_matrix, cancer_kmers["kmer"] == normal_matrix["kmer"], how='left_anti')
    
            logging.info("partitions: {}".format(cancer_kmers.rdd.getNumPartitions()))
            
            extension = '.tsv'
            path_final_fil = os.path.join(arg.output_dir, os.path.basename(cancer_path).split('.')[
                0] + '_ctlim{}_filt-normals-ctlim{}-{}samples'.format(arg.expr_limit_cancer, arg.expr_limit_normal, arg.n_samples_lim_normal) + extension)
            save_spark(cancer_kmers, arg.output_dir, path_final_fil, outpartitions=arg.out_partitions)


            ### Remove Uniprot
            logging.info("Filtering kmers in uniprot")
            cancer_kmers = remove_uniprot(spark, cancer_kmers, arg.uniprot, index_name)
            path_final_fil = os.path.join(arg.output_dir, os.path.basename(cancer_path).split('.')[
                0]  + '_ctlim{}_filt-normals-ctlim{}-{}_samples_filt-uniprot'.format(arg.expr_limit_cancer, arg.expr_limit_normal, arg.n_samples_lim_normal) + extension)
            save_spark(cancer_kmers, arg.output_dir, path_final_fil, outpartitions=arg.out_partitions)

            if os.path.exists(cancer_path_tmp):
                os.remove(cancer_path_tmp)


            #TODO Implement the intersection of the modelling tissues
            #TODO transfer to junction xpression

