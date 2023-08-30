import logging
import numpy as np
import pandas as pd
import sys
from pyspark.sql import functions as sf
from immunopepper.spark import split_only_found_annotation_backbone
from immunopepper.spark import filter_on_junction_kmer_annotated_flag


def process_libsize(path_lib, custom_normalizer):
    '''
    Loads and returns the normalisation values per sample
    Normalisation formulae: (count / 75 quantile expression in sample) * A
    If no normalisation factor is provided: A = median across samples
    If normalisation factor is provided: A = normalisation factor
    :param path_lib: str path with library size file
    :param custom_normalizer: custum normalisation factor
    :return: dataframe with 75 quantile expression in sample / A values
    '''
    lib = pd.read_csv(path_lib, sep='\t')
    if custom_normalizer:
        lib['libsize_75percent'] = lib['libsize_75percent'] / custom_normalizer
    else:
        lib['libsize_75percent'] = lib['libsize_75percent'] / np.median(lib['libsize_75percent'])
    lib['sample'] = [sample.replace('-', '').replace('.', '').replace('_','') for sample in lib['sample']]
    lib = lib.set_index('sample')
    return lib


def process_build_outputs(spark, index_name, coord_name, jct_col, jct_annot_col, rf_annot_col,
                        path_matrix=None, whitelist=None, cross_junction=None,
                        filterNeojuncCoord=None, filterAnnotatedRF=None,
                        output_dir=None, separate_back_annot=None, tot_batches=None, batch_id=None):
    '''
    Various preprocessing steps including
    - Loading of data from immunopepper mode build
    - Separation of the kmers derived solely from the annotation (expression is null in samples)
    - Filtering on junction and reading frame annotated status
    - Select whitelist samples
    :param spark: spark context
    :param index_name: str kmer column name
    :param jct_col: str junction column name
    :param jct_annot_col: str junction annotated flag column name
    :param rf_annot_col: str reading frame annotated flag column name
    :param path_matrix: str path for multi-sample count matrix from immunopepper mode build:  kmer|is_cross_junction|junctionAnnotated|readFrameAnnotated|sample_1|sample_2|...
    :param path_kmer_file: str path for single sample kmer file from immunopepper mode build:  kmer|id|segmentexpr|iscrossjunction|junctionexpr|junctionAnnotated|readFrameAnnotated|
    :param col_expr_kmer_file: str name of the expression column of interest in path_kmer_file. e.g. 'segmentexpr' or 'junctionexpr'
    :param target_sample: str sample for which the path_kmer_file has been computed
    :param whitelist: list whitelist for samples
    :param cross_junction: bool whether the count matrix contains junction counts or segment counts
    :param annot_flag: list with the intruction codes on how to treat the reading frame and junction annotated flags
    :param output_dir: str path for output directory
    :param filterNeojuncCoord: bool if True, filter for kmers from neojunctions,
     i.e. with non-annotated junction coordinates
    :param filterAnnotatedRF: bool if True, filter for kmers from annotated reading frames,
     i.e. with reading frames found in transcripts from the annotation, not propagated
    :param separate_back_annot: str, if provided, the kmers derived from the annotation only (i.e. without read support
    in any cohort sample) are saved into this path.
    - The kmers are thereby skipped from expression threshold filtering in the background matrix.
    - They will be by default removed from the foreground matrix as "annotated kmers"
    :param tot_batches: int batch mode only, total number of batches
    :param batch_id: int batch mode only, id of the batch
    :return: Preprocessed spark dataframe
    '''

    def RenamedCols(df):
        old_name = df.columns
        new_names = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in old_name]
        df = df.toDF(*new_names)
        return df


    def filter_whitelist(matrix, name_list, index_name, coord_name, jct_annot_col, rf_annot_col):
        name_list.extend([index_name, coord_name, jct_annot_col, rf_annot_col])
        return matrix.select([sf.col(name_) for name_ in name_list])


    if (path_matrix is not None):
        # Load immunopepper kmer candidates
        rename = False  # For development
        logging.info(f'Load input {path_matrix}')
        matrix = loader(spark, path_matrix, header=True)
        partitions_ = matrix.rdd.getNumPartitions()
        logging.info(f'...partitions: {partitions_}')

        if rename:
            logging.info("Rename")
            matrix = RenamedCols(matrix)

        # In batch mode: Dev remove kmers for the matrix based on the hash function
        if tot_batches:
            logging.info((f'Filter foreground and background based on Hash ; '
                          f'making {tot_batches} batches, select batch number {batch_id}'))

            matrix = matrix.filter(sf.abs(sf.hash('kmer') % tot_batches) == batch_id)

        # Filter on junction status depending on the content on the matrix
        matrix = matrix.drop(jct_col)

        # Separate kmers only present in the backbone annotation from the ones supported by the reads of any sample
        partitions_ = matrix.rdd.getNumPartitions()
        logging.info(f'...partitions: {partitions_}')
        if separate_back_annot:
            logging.info("Isolating kmers only in backbone annotation")
            matrix = split_only_found_annotation_backbone(separate_back_annot, output_dir, matrix, index_name,
                                                          jct_annot_col, rf_annot_col)

        # Filter according to annotation flag
        matrix = filter_on_junction_kmer_annotated_flag(matrix, jct_annot_col, rf_annot_col,
                                                      filterNeojuncCoord, filterAnnotatedRF)

        # Reduce samples (columns) to whitelist
        logging.info("Cast types")
        if whitelist is not None:
            whitelist = pd.read_csv(whitelist, sep='\t', header=None)[0].to_list()
            whitelist = [name_.replace('-', '').replace('.', '').replace('_', '') for name_ in whitelist]
            matrix = filter_whitelist(matrix, whitelist, index_name, coord_name, jct_annot_col, rf_annot_col)

        partitions_ = matrix.rdd.getNumPartitions()
        logging.info(f'...partitions: {partitions_}')
        return matrix
    else:
        return None


def remove_external_kmer_list(spark, path_external_kmer_database, spark_matrix, index_name, header=False):
    '''
    Takes an expression matrix of the format [kmers] x [samples, metadata] and removes an external database
    :param spark: spark context
    :param path_normal_kmer_list: str. path of file containing kmers to be removed
    :param spark_matrix: str path for count matrix
    :param index_name: str kmer column name
    :param header: bool. whether the file contains a header
    :return: Filtered matrix for external kmer database
    '''
    logging.info("Load input {}".format(path_external_kmer_database))
    external_database = loader(spark, path_external_kmer_database, header=header)
    external_database = external_database.select(sf.col(index_name))
    if spark_matrix:
        spark_matrix = spark_matrix.union(external_database)
    else:
        spark_matrix = external_database
    return spark_matrix


def loader(spark, path_kmer, header=False):
    '''
    Loads a parquet, csv or tsv (partitioned) file as a spark dataframe
    The user can provide either
    - a list of files with the same schema
    - a string. folder containing the files
    :param spark: spark context
    :param path_kmer: Either (1) a list of files or (b) a string with the basefolder where the files are placed
    :param header: bool. whether the file contains a header
    :return: loaded spark dataframe

    '''
    if ('parquet' in path_kmer) or ('pq' in path_kmer): #single file or folder
        matrix = spark.read.parquet(path_kmer)
    elif ('parquet' in path_kmer[0]) or ('pq' in path_kmer[0]): # list
        matrix = spark.read.parquet(*path_kmer)
    elif ('csv' in path_kmer[0]) or ('csv' in path_kmer):
        matrix = spark.read.csv(path_kmer, sep=',', header=header)
    elif ('tsv' in path_kmer[0]) or ('tsv' in path_kmer):
        matrix = spark.read.csv(path_kmer, sep=r'\t', header=header)
    elif ('gzip' in path_kmer[0]) or ('gzip' in path_kmer): #If tsv or csv not in path -> assume gzip are tab separated
        matrix = spark.read.csv(path_kmer, sep=r'\t', header=header)
    elif ('gz' in path_kmer[0]) or ('gz' in path_kmer):
        matrix = spark.read.csv(path_kmer, sep=r'\t', header=header)
    else:
        logging.error(f'Cannot determine file type of {path_kmer}. Please include .parquet .pq .tsv .csv suffix to the files or the folder (partitionned)')
        sys.exit()


    return matrix