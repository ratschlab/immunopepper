import logging
import os
from pyspark.sql import functions as sf


from .config import create_spark_session_from_config
from .config import default_spark_config
from .spark import combine_cancer
from .spark import combine_normals
from .spark import filter_expr_kmer
from .spark import filter_hard_threshold
from .spark import filter_statistical
from .spark import outlier_filtering
from .spark import pq_WithRenamedCols
from .spark import preprocess_kmer_file
from .spark import process_matrix_file
from .spark import process_libsize
from .spark import save_spark
from .spark import remove_uniprot



### Main
def mode_cancerspecif(arg):

    spark_cfg = default_spark_config(arg.cores, arg.mem_per_core, arg.parallelism)
    with create_spark_session_from_config(spark_cfg) as spark:
        # if os.path.exists(os.path.join(arg.output_dir, "checkpoint")):
        #     shutil.rmtree(os.path.join(arg.output_dir, "checkpoint"))
        # pathlib.Path(os.path.join(arg.output_dir, "checkpoint")).mkdir(exist_ok=True, parents=True)
        # spark.sparkContext.setCheckpointDir(os.path.join(arg.output_dir, "checkpoint"))

        index_name = 'kmer'
        jct_col = "iscrossjunction"
        extension = '.tsv'

        if arg.output_prefix:
            prefix =  arg.output_prefix + '_'
        else:
            prefix = ''

        if arg.paths_cancer_samples:
            if arg.expression_fields_c is None:
                expression_fields_orig = ['segment_expr', 'junction_expr']  # segment 1 , junction 2
            else:
                expression_fields_orig = arg.expression_fields_c
            expression_fields = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in
                                 expression_fields_orig]
            drop_cols = ['id']


        ### Preprocessing Libsize
        logging.info("\n \n >>>>>>>> Preprocessing libsizes")
        if arg.path_normal_libsize:
            libsize_n = process_libsize(arg.path_normal_libsize)
        else:
            libsize_n = None
        if arg.path_cancer_libsize:
            libsize_c = process_libsize(arg.path_cancer_libsize)
        else:
            libsize_c = None

        if arg.path_normal_kmer_list is None:
        ### Preprocessing Normals
            logging.info("\n \n >>>>>>>> Preprocessing Normal samples")

            normal_segm = process_matrix_file(spark, index_name, jct_col, arg.path_normal_matrix_segm, arg.output_dir, arg.whitelist_normal,
                            arg.parallelism, cross_junction=0)
            normal_junc = process_matrix_file(spark, index_name, jct_col, arg.path_normal_matrix_edge, arg.output_dir, arg.whitelist_normal,
                            arg.parallelism, cross_junction=1)
            normal_matrix = combine_normals(normal_segm, normal_junc, index_name)



            ### NORMALS: Statistical Filtering
            # Remove outlier kmers before statistical modelling (Very highly expressed kmers do not follow a NB, we classify them as expressed without hypothesis testing)
            if arg.statistical:
                logging.info("\n \n >>>>>>>> Normals: Perform statistical filtering")
                if arg.expr_high_limit_normal is not None:
                    logging.info( "... Normal kmers with expression >= {} in >= 1 sample are truly expressed. \n Will be substracted from cancer set".format(arg.expr_high_limit_normal))
                    high_expr_normals, normal_matrix = outlier_filtering(spark, normal_matrix, index_name, libsize_n, arg.expr_high_limit_normal)

                filter_statistical(spark, arg.tissue_grp_files, normal_matrix, index_name, arg.path_normal_matrix_segm,
                                       libsize_n, arg.threshold_noise, arg.output_dir, arg.cores)
            ### NORMALS: Hard Filtering
            else:
                logging.info("\n \n >>>>>>>> Normals: Perform Hard Filtering \n (expressed in {} samples with {} normalized counts)".format(arg.n_samples_lim_normal, arg.cohort_expr_support_normal))
                logging.info("expression filter")
                path_normal_kmers = filter_hard_threshold(normal_matrix, index_name, libsize_n, arg.output_dir, arg.cohort_expr_support_normal, arg.n_samples_lim_normal)
                normal_matrix = spark.read.csv(path_normal_kmers, sep=r'\t', header=False)
        else:
            logging.info("Load {}".format(arg.path_normal_kmer_list))
            if 'tsv' in arg.path_normal_kmer_list[0]:
                normal_matrix = spark.read.csv(arg.path_normal_kmer_list[0], sep=r'\t', header=False)
            else:
                normal_matrix = spark.read.parquet(*arg.path_normal_kmer_list)
        normal_matrix = normal_matrix.withColumnRenamed('_c0', index_name)
        normal_matrix = normal_matrix.select(sf.col(index_name))

        logging.info("\n \n >>>>>>>> Preprocessing Cancer samples")
        ### Apply filtering to foreground

        for cix, cancer_sample_ori in enumerate(arg.ids_cancer_samples):
            cancer_sample = cancer_sample_ori.replace('-', '').replace('.', '').replace('_', '')
            mutation_mode = arg.mut_cancer_samples[cix]

            ## Cancer file is matrix
            if arg.path_cancer_matrix_segm or arg.path_cancer_matrix_edge:
                mutation_mode = arg.mut_cancer_samples[0]
                # Preprocess cancer samples
                cancer_segm = process_matrix_file(spark, index_name, jct_col, arg.path_cancer_matrix_segm,
                                                  arg.output_dir, arg.whitelist_cancer,
                                                  arg.parallelism, cross_junction=0)
                cancer_junc = process_matrix_file(spark, index_name, jct_col, arg.path_cancer_matrix_edge,
                                                  arg.output_dir, arg.whitelist_cancer,
                                                  arg.parallelism, cross_junction=1)
                cancer_matrix = combine_cancer(cancer_segm, cancer_junc, index_name)
                # Apply expression filter to foreground
                if arg.scratch_dir:
                    cancer_out = os.environ[arg.scratch_dir]
                else:
                    cancer_out = arg.output_dir
                #cross sample filter
                path_cancer_kmers = filter_hard_threshold(cancer_matrix, index_name, libsize_c, cancer_out,
                                                          arg.cohort_expr_support_cancer, arg.n_samples_lim_cancer,
                                                          target_sample=cancer_sample,
                                                          tag='cancer_{}'.format(mutation_mode))
                valid_foreground = spark.read.csv(path_cancer_kmers, sep=r'\t', header=False)
                valid_foreground = valid_foreground.withColumnRenamed('_c0', index_name)
                valid_foreground = valid_foreground.select(sf.col(index_name))
                cancer_cross_filter = cancer_matrix.join(valid_foreground, ["kmer"], # Probably do union differently
                                                   how='right').select([index_name, cancer_sample])
                # sample specific filter
                cancer_sample_filter = cancer_matrix.select([index_name, cancer_sample])
                cancer_sample_filter = filter_expr_kmer(cancer_sample_filter, cancer_sample, arg.sample_expr_support_cancer)
                
                cancer_kmers = cancer_cross_filter.union(cancer_sample_filter).distinct()
            ## Cancer file is kmer file
            if arg.paths_cancer_samples:
                cancer_path = arg.paths_cancer_samples[cix]
                rename = True # development
                if rename: 
                    cancer_path_tmp = pq_WithRenamedCols(cancer_path, arg.output_dir)
                else: 
                    cancer_path_tmp = cancer_path
                cancer_kmers = spark.read.parquet(cancer_path_tmp)
                # Preprocess cancer samples
                cancer_junc = preprocess_kmer_file(cancer_kmers, cancer_sample, drop_cols,expression_fields, jct_col,
                                                                                        index_name, libsize_c, 1)
                cancer_segm = preprocess_kmer_file(cancer_kmers, cancer_sample, drop_cols, expression_fields, jct_col,
                                                                                        index_name, libsize_c, 0)
                # Apply expression filter to foreground
                cancer_junc = filter_expr_kmer(cancer_junc, expression_fields[1], arg.sample_expr_support_cancer)
                cancer_segm = filter_expr_kmer(cancer_segm, expression_fields[0], arg.sample_expr_support_cancer)
                cancer_kmers = combine_cancer(cancer_junc, cancer_segm, index_name)


        ## Cancer \ normals
            logging.info("\n \n >>>>>>>> Cancers: Perform differential filtering sample {}".format(cancer_sample_ori))
            logging.info("partitions: {}".format(cancer_kmers.rdd.getNumPartitions()))

            # outpaths
            base_path_final= os.path.join(arg.output_dir, '{}{}_{}_SampleLim{}CohortLim{}Across{}_FiltNormalsCohortlim{}Across{}'.format(prefix,
                                                           cancer_sample_ori, mutation_mode, arg.sample_expr_support_cancer,
                                                                        arg.cohort_expr_support_cancer, arg.n_samples_lim_cancer,
                                                                        arg.cohort_expr_support_normal, arg.n_samples_lim_normal))
            path_filter_final = base_path_final + extension
            path_filter_final_uniprot  = base_path_final + '_FiltUniprot'+ extension

            # Remove background from foreground
            logging.info("Filtering normal background")
            cancer_kmers = cancer_kmers.join(normal_matrix, cancer_kmers["kmer"] == normal_matrix["kmer"], how='left_anti')
            logging.info("partitions: {}".format(cancer_kmers.rdd.getNumPartitions()))
            save_spark(cancer_kmers, arg.output_dir, path_filter_final, outpartitions=arg.out_partitions)

            # Remove Uniprot
            logging.info("Filtering kmers in uniprot")
            cancer_kmers = remove_uniprot(spark, cancer_kmers, arg.uniprot, index_name)
            save_spark(cancer_kmers, arg.output_dir, path_filter_final_uniprot, outpartitions=arg.out_partitions)


            if arg.paths_cancer_samples and os.path.exists(cancer_path_tmp) and rename:
                os.remove(cancer_path_tmp)


            #TODO Implement the intersection of the modelling tissues

