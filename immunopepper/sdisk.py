import logging
import numpy as np
import os
import pathlib

def save_spark(cancer_kmers, output_dir, path_final_fil, outpartitions=None):
    '''
    Saves a spark dataframe as a single or partitioned csv file
    :param cancer_kmers: spark dataframe matrix with expression counts for cancer
    :param output_dir: str path for output directory
    :param path_final_fil: str path to save the spark dataframe
    :param outpartitions: int number of partitions for saving
    '''
    # save
    logging.info(f'>>>> Save to {path_final_fil}')
    pathlib.Path(output_dir).mkdir(exist_ok=True, parents=True)
    if outpartitions is not None:
        cancer_kmers.repartition(outpartitions).write.mode('overwrite')\
            .options(header="true", sep="\t", compression="gzip").format("tsv.gz").csv(path_final_fil)
    else:
        cancer_kmers.write.mode('overwrite')\
            .options(header="true", sep="\t", compression="gzip").format("tsv.gz").csv(path_final_fil)


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
        header = (f'{"sample"}\t{"mutation_mode"}\t{"min_sample_reads"}\t{"#_of_cohort_samples"}\t'
                  f'{"reads_per_cohort_sample"}\t{"#_normal_samples_allowed"}\t{"normal_cohort_id"}'
                  f'\t{"reads_per_normal_sample"}')
        line =   (f'{cancer_sample_ori}\t{mutation_mode}\t{sample_expr_support_cancer}\t{n_samples_lim_cancer}'
                  f'\t{cohort_expr_support_cancer}\t{n_samples_lim_normal}\t{id_normals}'
                  f'\t{cohort_expr_support_normal}')

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


def redirect_interm(interm_dir_norm, interm_dir_canc, output_dir):
    '''
    Set the directory to save intermediary file
    - The output directory
    - Any other specified normal or cancer directory
    Default. Uses output directory
    :param interm_dir_norm: str custom scatch dir path to save intermediate normal files
    :param interm_dir_canc: str custom scatch dir path to save intermediate cancer files
    :param output_dir: str output directory for the filtered matrix
    :return:
    '''

    if interm_dir_canc:
        cancer_out = interm_dir_canc
    else:
        cancer_out = output_dir
    if interm_dir_norm:
        normal_out = interm_dir_norm
    else:
        normal_out = output_dir
    return normal_out, cancer_out


def check_interm_files(out_dir, expr_limit, n_samples_lim, target_sample='', tag='normals', batch_tag=''):
    '''
    Filtering steps for normal (resp. cancer) samples are saved as intermediate files because it is an expensive operation
    The function checks the presence of the intermediate filtering files to decide whether to perform
    - the full preprocessing + threshold filtering steps
    - or simply re-load the intermediate files
    :param out_dir: str path for output directory
    :param expr_limit: float expression limit threshold to keep a kmer
    :param n_samples_lim: int number of samples that need to pass the expression limit
    :param target_sample: str name of the sample of interest.
    To be excluded in the number of samples that pass the expression limit
    :param tag: str tag related to the type of samples. Example cancer or normal
    :param batch_tag: str batch mode, batch tag to be appended to intermediate file
    :returns:
    - launch_preprocess: bool, whether to perform the full preprocessing + threshold filtering steps
     or simply re-load the intermediate files
    - path_interm_matrix_for_express_threshold, path_interm_matrix_for_sample_threshold,
     path_interm_kmers_annotOnly are respectively the path (str) where
     the expression-filtered matrix, the sample-filtered matrix and
     the kmers derived solely from the annotation are saved
    '''
    base_n_samples = 1
    base_expr = 0.0
    format_tag = '.tsv.gz'
    # For cancer matrix intermediate file the recurrence filter is not applied to the target sample
    if target_sample:
        suffix = f'Except{target_sample}'
    else:
        suffix = ''

    # For normal samples the expression threshold filtering is not applied to the kmers found only in the annotation
    # but not in the background samples. These kmers will be by default removed from the foreground matrix.
    if tag == 'normals':
        path_interm_kmers_annotOnly = os.path.join(out_dir, f'kmers_derived_solely_from_annotation{format_tag}')
    else:
        path_interm_kmers_annotOnly = None

    # Saving paths

    path_interm_matrix_for_express_threshold = os.path.join(out_dir,
                          f'interm_{tag}_combiExprCohortLim{expr_limit}Across{base_n_samples}{suffix}{batch_tag}{format_tag}')
    path_interm_matrix_for_sample_threshold = os.path.join(out_dir,
                              f'interm_{tag}_combiExprCohortLim{base_expr}Across{base_n_samples}{suffix}{batch_tag}{format_tag}')
    # Check existence
    if (expr_limit and os.path.isfile(os.path.join(path_interm_matrix_for_express_threshold, '_SUCCESS'))) \
        and (n_samples_lim is not None and os.path.isfile(os.path.join(path_interm_matrix_for_sample_threshold, '_SUCCESS'))):

        logging.info((f'Intermediate {tag} filtering already performed in: {path_interm_matrix_for_express_threshold} '
                      f' and {path_interm_matrix_for_sample_threshold}. Re-loading {tag} intermediate data...'))
        logging.info((f'Proceed with care! Using intermediate files means ignoring --filterNeojuncCoord, '
                      f'--filterAnnotatedRF parameter.'))
        launch_preprocess = False
    else:
        logging.info(f'At least one intermediate {tag} filtering file is missing.')
        logging.info(f'Will compute full filtering steps according to user input parameters')
        launch_preprocess = True

    return launch_preprocess, \
           path_interm_matrix_for_express_threshold, \
           path_interm_matrix_for_sample_threshold, \
           path_interm_kmers_annotOnly



def filtered_path(arg, cancer_sample_ori, mutation_mode, normal_files, batch_tag, extension):
    '''
    :param arg: argument class from argparse
    :param cancer_sample_ori: str. name of the cancer sample
    :param mutation_mode: str. mutation mode flag
    :param normal_files: bool. whether the normal files are used as an input
    :param batch_tag: str. tag for the batch
    :param extension: str. saving format extension
    :return:
    path_filter_final
    path_filter_final_uniprot
    '''
    base_path_final = os.path.join(arg.output_dir,
                                   (f'{arg.tag_prefix}{cancer_sample_ori}_{mutation_mode}_'
                                    f'SampleLim{arg.sample_expr_support_cancer}'
                                    f'CohortLim{arg.cohort_expr_support_cancer}'
                                    f'Across{arg.n_samples_lim_cancer}'))
    if normal_files:
        base_path_final += (f'_FiltNormals{arg.tag_normals}'
                            f'Cohortlim{arg.cohort_expr_support_normal}'
                            f'Across{arg.n_samples_lim_normal}')

    path_filter_final = base_path_final + batch_tag + extension
    path_filter_final_uniprot = base_path_final + '_FiltUniprot' + batch_tag + extension

    return path_filter_final, path_filter_final_uniprot