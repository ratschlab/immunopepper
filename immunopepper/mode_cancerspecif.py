import logging
import os

from immunopepper.sdisk import check_interm_files
from immunopepper.sdisk import filtered_path
from immunopepper.sdisk import output_count
from immunopepper.sdisk import redirect_interm
from immunopepper.sdisk import save_output_count
from immunopepper.sdisk import save_spark
from immunopepper.sloaders import apply_preprocess
from immunopepper.sloaders import inputs_to_modes
from immunopepper.sloaders import process_libsize
from immunopepper.sloaders import remove_external_kmer_list
from immunopepper.spark import combine_hard_threshold_cancers
from immunopepper.spark import combine_hard_threshold_normals
from immunopepper.spark import filter_expr_kmer
from immunopepper.spark import filter_hard_threshold
from immunopepper.spark import remove_uniprot
from immunopepper.spark_config import create_spark_session_from_config
from immunopepper.spark_config import default_spark_config
from pyspark.sql.functions import broadcast




### Main
def mode_cancerspecif(arg):

    spark_cfg = default_spark_config(arg.cores, arg.mem_per_core, arg.parallelism, tmp_dir=arg.scratch_dir)
    with create_spark_session_from_config(spark_cfg) as spark:
        # if os.path.exists(os.path.join(arg.output_dir, "checkpoint")):
        #     shutil.rmtree(os.path.join(arg.output_dir, "checkpoint"))
        # pathlib.Path(os.path.join(arg.output_dir, "checkpoint")).mkdir(exist_ok=True, parents=True)
        # spark.sparkContext.setCheckpointDir(os.path.join(arg.output_dir, "checkpoint"))

        ### Some variable definitions
        index_name = 'kmer'
        coord_name = 'coord'
        jct_col = "iscrossjunction"
        jct_annot_col = "junctionAnnotated"
        rf_annot_col = "readFrameAnnotated"
        extension = '.tsv.gz'
        if arg.batch_id is not None:
            batch_tag = '_batch{}_{}'.format(arg.batch_id, arg.tot_batches)
        else:
            batch_tag = ''
        if (arg.tag_prefix) and arg.tag_prefix[-1] != '_':
            arg.tag_prefix = arg.tag_prefix + '_'
        report_count = []
        report_steps = []
        normal_out, cancer_out = redirect_interm(arg.interm_dir_norm,
                                                  arg.interm_dir_canc, arg.output_dir)
        recurrence_cancer, recurrence_normal, cancer_files, normal_files = inputs_to_modes(arg)

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

        launch_preprocess_normal, path_normal_for_express_threshold, path_normal_for_sample_threshold, \
        path_interm_kmers_annotOnly = check_interm_files(normal_out, arg.cohort_expr_support_normal,
                                                         arg.n_samples_lim_normal, tag='normals', batch_tag=batch_tag)
        ### Preprocessing Normals
        if normal_files:
            if launch_preprocess_normal: # else do not need to launch because intermediate files are present
                logging.info("\n \n >>>>>>>> Preprocessing Normal samples")

                normal_matrix = apply_preprocess(spark, 'N', index_name, coord_name, jct_col, jct_annot_col, rf_annot_col,
                        arg.whitelist_normal, arg.path_normal_matrix_segm, arg.path_normal_matrix_edge,
                        arg.filterNeojuncCoord, arg.filterAnnotatedRF, arg.tot_batches, arg.batch_id,
                        normal_out, path_interm_kmers_annotOnly)


            # Hard Filtering
            if recurrence_normal:
                logging.info((f'\n \n >>>>>>>> Normals: Perform Hard Filtering \n '
                              f'(expressed in {arg.n_samples_lim_normal} samples'
                              f' with {arg.cohort_expr_support_normal} normalized counts'))
                logging.info("expression filter")
                inter_matrix_expr, inter_matrix_sample = filter_hard_threshold(normal_matrix, index_name, coord_name, jct_annot_col,
                                                                               rf_annot_col, libsize_n,
                                                                               arg.cohort_expr_support_normal,
                                                                               arg.n_samples_lim_normal,
                                                                               path_normal_for_express_threshold,
                                                                               path_normal_for_sample_threshold,
                                                                               on_the_fly=arg.on_the_fly, tag='normals')

                normal_matrix = combine_hard_threshold_normals(spark, path_normal_for_express_threshold,
                                                               path_normal_for_sample_threshold,
                                                               inter_matrix_expr, inter_matrix_sample,
                                                               arg.n_samples_lim_normal, index_name)

                # Add back kmer annot
                normal_matrix = remove_external_kmer_list(spark, path_interm_kmers_annotOnly,
                                                      normal_matrix, index_name, header=True)

        # Additional kmer backgrounds filtering
        if arg.path_normal_kmer_list is not None:
            normal_matrix = remove_external_kmer_list(spark, arg.path_normal_kmer_list, normal_matrix, index_name, header=True)

        ### Apply filtering to foreground
        for cix, cancer_sample_ori in enumerate(arg.ids_cancer_samples):
            cancer_sample = cancer_sample_ori.replace('-', '').replace('.', '').replace('_', '')
            mutation_mode = arg.mut_cancer_samples[cix]

            launch_filter_cancer, path_cancer_for_express_threshold, path_cancer_for_sample_threshold, _ \
                = check_interm_files(cancer_out, arg.cohort_expr_support_cancer, arg.n_samples_lim_cancer, \
                                     target_sample=cancer_sample, tag=f'cancer_{mutation_mode}', batch_tag=batch_tag)

            ## Cancer file Filters
            if cancer_files:
                logging.info("\n \n >>>>>>>> Preprocessing Cancer sample {}  ".format(cancer_sample_ori))
                mutation_mode = arg.mut_cancer_samples[0]

                cancer_matrix = apply_preprocess(spark, 'C', index_name, coord_name, jct_col, jct_annot_col, rf_annot_col,
                        arg.whitelist_cancer, arg.path_cancer_matrix_segm, arg.path_cancer_matrix_edge,
                        arg.filterNeojuncCoord, arg.filterAnnotatedRF, arg.tot_batches, arg.batch_id)

                # cancer sample-specific filter

                cancer_sample_filter = cancer_matrix.select([index_name, coord_name, cancer_sample, jct_annot_col, rf_annot_col])

                cancer_sample_filter = filter_expr_kmer(cancer_sample_filter, cancer_sample, 0) #Keep kmers expressed
                # Counting step: Retrieve initial number of kmers in sample
                output_count(arg.output_count, cancer_sample_filter, report_count, report_steps, 'Init_cancer')

                if arg.output_count and (arg.sample_expr_support_cancer != 0):
                    cancer_sample_filter = filter_expr_kmer(cancer_sample_filter, cancer_sample,
                                                            arg.sample_expr_support_cancer, libsize_c) #Keep kmers expressed >= threshold

                    # Counting step: Retrieve number of kmers in sample after filtering on cancer expression
                    output_count(arg.output_count, cancer_sample_filter, report_count, report_steps, 'Filter_Sample')

                else:
                    output_count(arg.output_count, cancer_sample_filter, report_count, report_steps, 'Filter_Sample')

                # cancer cross-cohort filter
                if recurrence_cancer:
                    inter_matrix_expr_c, inter_matrix_sample_c = filter_hard_threshold(cancer_matrix, index_name,coord_name,
                                                                                       jct_annot_col, rf_annot_col,
                                                                                       libsize_c,
                                                                                       arg.cohort_expr_support_cancer,
                                                                                       arg.n_samples_lim_cancer,
                                                                                       path_cancer_for_express_threshold,
                                                                                       path_cancer_for_sample_threshold,
                                                                                       target_sample=cancer_sample,
                                                                                       on_the_fly=arg.on_the_fly,
                                                                                       tag=f'cancer_{mutation_mode}')
                    cancer_cross_filter = combine_hard_threshold_cancers(spark, cancer_matrix,
                                                                         path_cancer_for_express_threshold,
                                                                         path_cancer_for_sample_threshold,
                                                                         inter_matrix_expr_c, inter_matrix_sample_c,
                                                                         arg.cohort_expr_support_cancer,
                                                                         arg.n_samples_lim_cancer,
                                                                         index_name)

                    cancer_cross_filter = cancer_cross_filter.select([index_name, coord_name, cancer_sample,
                                                                      jct_annot_col, rf_annot_col])
                    if arg.cancer_support_union:
                        logging.info("support union")
                        cancer_kmers = cancer_cross_filter.union(cancer_sample_filter).distinct()
                    else:
                        logging.info("support intersect")
                        cancer_kmers = cancer_cross_filter.join(cancer_sample_filter.select([index_name]), ["kmer"],
                                                                how='inner')
                else:
                    cancer_kmers = cancer_sample_filter


                output_count(arg.output_count, cancer_kmers, report_count, report_steps, 'Filter_Sample_Cohort')

            # Outpaths
            path_filter_final, path_filter_final_uniprot = filtered_path(arg, cancer_sample_ori, mutation_mode,
                                                                         recurrence_normal, batch_tag, extension)
            # Remove Background
            if normal_files:
                logging.info("\n \n >>>>>>>> Cancers: Perform differential filtering")
                partitions_ = cancer_kmers.rdd.getNumPartitions()
                logging.info(f'partitions: {partitions_}')
                logging.info("Filtering normal background")
                cancer_kmers = broadcast(cancer_kmers)
                cancer_kmers = cancer_kmers.join(normal_matrix, cancer_kmers["kmer"] == normal_matrix["kmer"],
                                             how='left_anti')

            # Save main output
            save_spark(cancer_kmers, arg.output_dir, path_filter_final, outpartitions=arg.out_partitions)
            output_count(arg.output_count, cancer_kmers, report_count, report_steps,
                         'Filter_Sample_Cohort_CohortNormal')

            # Remove Uniprot #TODO redondant with external database
            if arg.uniprot is not None:
                logging.info("Filtering kmers in uniprot")
                cancer_kmers = remove_uniprot(spark, cancer_kmers, arg.uniprot, index_name)
                save_spark(cancer_kmers, arg.output_dir, path_filter_final_uniprot, outpartitions=arg.out_partitions)
                output_count(arg.output_count, cancer_kmers, report_count, report_steps,
                         'Filter_Sample_Cohort_CohortNormal_Uniprot')

            save_output_count(arg.output_count, report_count, report_steps, arg.tag_normals,
                                cancer_sample_ori, mutation_mode, arg.sample_expr_support_cancer,
                                arg.cohort_expr_support_cancer, arg.n_samples_lim_cancer,
                                arg.cohort_expr_support_normal, arg.n_samples_lim_normal, arg.tag_normals)




