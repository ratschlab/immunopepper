import logging
import os
from pyspark.sql import functions as sf


from .config import create_spark_session_from_config
from .config import default_spark_config
from .spark import combine_cancer
from .spark import combine_hard_threshold_cancers
from .spark import combine_hard_threshold_normals
from .spark import combine_normals
from .spark import filter_expr_kmer
from .spark import filter_hard_threshold
from .spark import filter_statistical
from .spark import loader
from .spark import outlier_filtering
from .spark import output_count
from .spark import pq_WithRenamedCols
from .spark import preprocess_kmer_file
from .spark import process_matrix_file
from .spark import process_libsize
from .spark import redirect_scratch
from .spark import remove_uniprot
from .spark import save_output_count
from .spark import save_spark





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
        if arg.batch_id is not None:
            batch_tag = '_batch{}'.format(arg.batch_id)
        else:
            batch_tag = ''
        if (arg.tag_prefix is not None) and arg.tag_prefix[-1] != '_':
            arg.tag_prefix = arg.tag_prefix + '_'
        report_count = []
        report_steps = []

        normal_out, cancer_out = redirect_scratch(arg.scratch_dir, arg.interm_dir_norm, arg.interm_dir_canc, arg.output_dir)


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
            libsize_n = process_libsize(arg.path_normal_libsize, arg.normalizer_normal_libsize)
        else:
            libsize_n = None
        if arg.path_cancer_libsize:
            libsize_c = process_libsize(arg.path_cancer_libsize, arg.normalizer_cancer_libsize)
        else:
            libsize_c = None

        normal_matrix = None
        ### Preprocessing Normals
        if (arg.path_normal_matrix_segm is not None) or (arg.path_normal_matrix_edge is not None):
            logging.info("\n \n >>>>>>>> Preprocessing Normal samples")

            normal_segm = process_matrix_file(spark, index_name, jct_col, arg.path_normal_matrix_segm, arg.output_dir, arg.whitelist_normal,
                            arg.parallelism, cross_junction=0, tot_batches=arg.tot_batches, batch_id=arg.batch_id)
            normal_junc = process_matrix_file(spark, index_name, jct_col, arg.path_normal_matrix_edge, arg.output_dir, arg.whitelist_normal,
                            arg.parallelism, cross_junction=1, tot_batches=arg.tot_batches, batch_id=arg.batch_id)
            normal_matrix = combine_normals(normal_segm, normal_junc, index_name)

            # NORMALS: Statistical Filtering # Remove outlier kmers before statistical modelling (Very highly expressed kmers do not follow a NB, we classify them as expressed without hypothesis testing)
            if arg.statistical:
                logging.info("\n \n >>>>>>>> Normals: Perform statistical filtering")
                if arg.expr_high_limit_normal is not None:
                    logging.info( "... Normal kmers with expression >= {} in >= 1 sample are truly expressed. \n Will be substracted from cancer set".format(arg.expr_high_limit_normal))
                    high_expr_normals, normal_matrix = outlier_filtering(spark, normal_matrix, index_name, libsize_n, arg.expr_high_limit_normal)

                filter_statistical(spark, arg.tissue_grp_files, normal_matrix, index_name, arg.path_normal_matrix_segm,
                                       libsize_n, arg.threshold_noise, arg.output_dir, arg.cores)
            # NORMALS: Hard Filtering
            else:
                logging.info("\n \n >>>>>>>> Normals: Perform Hard Filtering \n (expressed in {} samples with {} normalized counts)".format(arg.n_samples_lim_normal, arg.cohort_expr_support_normal))
                logging.info("expression filter")
                path_normal_kmers_e, path_normal_kmers_s = filter_hard_threshold(normal_matrix, index_name, libsize_n, normal_out, arg.cohort_expr_support_normal, arg.n_samples_lim_normal, batch_tag=batch_tag)
                normal_matrix = combine_hard_threshold_normals(spark, path_normal_kmers_e, path_normal_kmers_s, arg.n_samples_lim_normal, index_name)

        if arg.path_normal_kmer_list is not None:
            logging.info("Load {}".format(arg.path_normal_kmer_list))
            normal_list = loader(spark, arg.path_normal_kmer_list)
            normal_list = normal_list.select(sf.col(index_name))
            if normal_matrix:
                normal_matrix = normal_matrix.union(normal_list).distinct()
            else:
                normal_matrix = normal_list


        ### Apply filtering to foreground

        for cix, cancer_sample_ori in enumerate(arg.ids_cancer_samples):
            cancer_sample = cancer_sample_ori.replace('-', '').replace('.', '').replace('_', '')
            mutation_mode = arg.mut_cancer_samples[cix]
            logging.info("\n \n >>>>>>>> Preprocessing Cancer sample {}  ".format(cancer_sample_ori))

            ## Cancer file is matrix
            if arg.path_cancer_matrix_segm or arg.path_cancer_matrix_edge:
                mutation_mode = arg.mut_cancer_samples[0]
                # Preprocess cancer samples
                cancer_segm = process_matrix_file(spark, index_name, jct_col, arg.path_cancer_matrix_segm,
                                                  arg.output_dir, arg.whitelist_cancer,
                                                  arg.parallelism, cross_junction=0, tot_batches=arg.tot_batches, batch_id=arg.batch_id)
                cancer_junc = process_matrix_file(spark, index_name, jct_col, arg.path_cancer_matrix_edge,
                                                  arg.output_dir, arg.whitelist_cancer,
                                                  arg.parallelism, cross_junction=1, tot_batches=arg.tot_batches, batch_id=arg.batch_id)
                cancer_matrix = combine_cancer(cancer_segm, cancer_junc, index_name)


                # sample specific filter
                cancer_sample_filter = cancer_matrix.select([index_name, cancer_sample])
                cancer_sample_filter = filter_expr_kmer(cancer_sample_filter, cancer_sample, 0) #Retrieve initial number of kmers in sample
                output_count(arg.output_count, cancer_sample_filter, report_count, report_steps, 'Init_cancer')

                if arg.output_count and (arg.sample_expr_support_cancer != 0):
                    cancer_sample_filter = cancer_matrix.select([index_name, cancer_sample])
                    cancer_sample_filter = filter_expr_kmer(cancer_sample_filter, cancer_sample,
                                                            arg.sample_expr_support_cancer)
                    output_count(arg.output_count, cancer_sample_filter, report_count, report_steps, 'Filter_Sample')

                else:
                    output_count(arg.output_count, cancer_sample_filter, report_count, report_steps, 'Filter_Sample')


                #cross sample filter
                if (arg.cohort_expr_support_cancer is not None) and (arg.n_samples_lim_cancer is not None):
                    path_cancer_kmers_e, path_cancer_kmers_s  = filter_hard_threshold(cancer_matrix, index_name, libsize_c, cancer_out,
                                                              arg.cohort_expr_support_cancer, arg.n_samples_lim_cancer,
                                                              target_sample=cancer_sample,
                                                              tag='cancer_{}'.format(mutation_mode),
                                                              batch_tag=batch_tag)
                    cancer_cross_filter = combine_hard_threshold_cancers(spark, cancer_matrix, path_cancer_kmers_e, arg.cohort_expr_support_cancer, arg.n_samples_lim_cancer, index_name, cancer_sample)


                    if arg.cancer_support_union:
                        cancer_kmers = cancer_cross_filter.union(cancer_sample_filter).distinct()
                    else:
                        cancer_kmers = cancer_cross_filter.join(cancer_sample_filter.drop(cancer_sample), ["kmer"], how='inner')
                else:
                    cancer_kmers = cancer_sample_filter


                output_count(arg.output_count, cancer_kmers, report_count, report_steps, 'Filter_Sample_Cohort')

            ## Cancer file is kmer file
            if arg.paths_cancer_samples:
                cancer_path = arg.paths_cancer_samples[cix]
                rename = True # development
                if rename:
                    cancer_kmers = pq_WithRenamedCols(spark, cancer_path, arg.output_dir)
                else:
                    cancer_kmers = spark.read.parquet(cancer_path)

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
            logging.info("\n \n >>>>>>>> Cancers: Perform differential filtering")
            logging.info("partitions: {}".format(cancer_kmers.rdd.getNumPartitions()))

            # outpaths
            base_path_final= os.path.join(arg.output_dir, '{}{}_{}_SampleLim{}CohortLim{}Across{}_FiltNormals{}Cohortlim{}Across{}'.format(
                                                            arg.tag_prefix, cancer_sample_ori, mutation_mode, arg.sample_expr_support_cancer,
                                                            arg.cohort_expr_support_cancer, arg.n_samples_lim_cancer,
                                                            arg.tag_normals, arg.cohort_expr_support_normal, arg.n_samples_lim_normal))

            path_filter_final = base_path_final + batch_tag + extension
            path_filter_final_uniprot  = base_path_final + '_FiltUniprot'+ batch_tag + extension

            # Remove background from foreground
            logging.info("Filtering normal background")
            cancer_kmers = cancer_kmers.join(normal_matrix, cancer_kmers["kmer"] == normal_matrix["kmer"], how='left_anti')
            logging.info("partitions: {}".format(cancer_kmers.rdd.getNumPartitions()))
            save_spark(cancer_kmers, arg.output_dir, path_filter_final, outpartitions=arg.out_partitions)
            output_count(arg.output_count, cancer_kmers, report_count, report_steps,
                         'Filter_Sample_Cohort_CohortNormal')

            # Remove Uniprot
            logging.info("Filtering kmers in uniprot")
            cancer_kmers = remove_uniprot(spark, cancer_kmers, arg.uniprot, index_name)
            save_spark(cancer_kmers, arg.output_dir, path_filter_final_uniprot, outpartitions=arg.out_partitions)
            output_count(arg.output_count, cancer_kmers, report_count, report_steps,
                         'Filter_Sample_Cohort_CohortNormal_Uniprot')

            save_output_count(arg.output_count, report_count, report_steps, arg.tag_normals,
                                cancer_sample_ori, mutation_mode, arg.sample_expr_support_cancer,
                                arg.cohort_expr_support_cancer, arg.n_samples_lim_cancer,
                                arg.cohort_expr_support_normal, arg.n_samples_lim_normal, arg.tag_normals)

            #if arg.paths_cancer_samples and os.path.exists(cancer_path_tmp) and rename:
            #    os.remove(cancer_path_tmp)


            #TODO Implement the intersection of the modelling tissues

