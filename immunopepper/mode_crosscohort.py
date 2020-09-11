from collections import defaultdict
import glob
import logging
import os
import pathlib
from pyspark.sql import functions as sf
from pyspark.sql import types as st
import shutil
import sys
from typing import Iterable

from .config import create_spark_session
from .config import create_spark_session_from_config
from .config import default_spark_config

def organize_inputs(input_dir, backgr_suffix, kmer, samples, mutation_modes, compres_suffix, skip_reorganize):
    to_combine = []
    combined_dir = "junction_kmers"
    for sample in samples:
        sample_basedir = os.path.join(input_dir, sample)
        pathlib.Path(os.path.join(sample_basedir, "junction_kmers")).mkdir(exist_ok= True, parents= True)
        jct_inputs = ['{}_junction_{}mer{}.pq{}'.format(mutation, kmer, backgr_suffix, compres_suffix)
                      for mutation in mutation_modes]
        if not skip_reorganize:
            for input_file in jct_inputs:
                source = os.path.join(sample_basedir, input_file)
                target = os.path.join(sample_basedir, combined_dir, input_file)
                if (not os.path.isfile(source)) and os.path.isfile(target):
                    logging.error("File {} seems to have been moved to {}, rerun with --skip-filegrouping".format(source, target))
                    sys.exit(1)
                elif os.path.isfile(source) and os.path.isfile(target):
                    logging.warning("Overwriting {} with {} file".format(target, source))
                    shutil.move(source, target)
                elif (not os.path.isfile(source)) and (not os.path.isfile(target)):
                    logging.error("Inputs {} not generated, run mode_build".format(source))
                    sys.exit(1)
                else:
                    logging.info("Moving {} to {}".format(source, target))
                    shutil.move(source, target)
        to_combine.append((os.path.join(sample_basedir, combined_dir), os.path.join(sample_basedir, "combined_junction_kmers.pq")))
    return to_combine


def run_combine(inf: Iterable[str],
                outf: str,
                partition_field: str,
                grp_fields: Iterable[str],
                expression_fields: Iterable[str],
                max_aggregation_fields: Iterable[str],
                cores: int, memory_per_core:int, partitions: int):
    spark_cfg = default_spark_config(cores, memory_per_core)

    with create_spark_session_from_config(spark_cfg) as spark:
        df = spark.read.parquet(inf + '/') #TODO check
        custom_max = sf.udf(collapse_values, st.StringType())
        for e_col in expression_fields:
            df = df.withColumn(e_col, custom_max(e_col))
        agg_cols = [sf.max(c).alias('max_' + c) for c in max_aggregation_fields]

        df_g = (df.repartition(partitions, partition_field).
                        groupBy(grp_fields).
                        agg(*agg_cols)
               )

        df_g.write.mode('overwrite').parquet(outf)


def run_pivot(input_, output_dir, output_suffix, sample_list, aggregation_col, partitions, cores, memory_per_core):
    sample_set = set(sample_list)
    for expression_col in aggregation_col:
        output_path = os.path.join(output_dir, "crosssamples_" + expression_col + output_suffix + '.pq')
        file_info = defaultdict(list, {})
        for pq_folder in input_:
            for partition in glob.glob(pq_folder + '/*'):
                index_partition = '-'.join(partition.split('/')[-1].split('-')[0:2])
                if index_partition != '_SUCCESS':
                    file_info[index_partition].append(partition)
        print(file_info)

        with create_spark_session(cores, memory_per_core) as spark:
            for part_id in file_info:
                print(f"working on part {part_id}")
                part_files = file_info[part_id]

                df = spark.read.parquet(*part_files)
                print(df.show())

                path_parts = sf.split(sf.input_file_name(), '/')
                df = df.withColumn('sample', path_parts.getItem(
                    sf.size(path_parts) - 3))  # depends on the sample position in the path
                print(df.show())

                df_junction = (df.groupBy('kmer').
                               agg(sf.max('max_is_cross_junction').alias('is_cross_junction')).
                               withColumnRenamed('kmer', 'kmer_junction')
                               )
                print(df_junction.show())

                print(expression_col)
                df_pivoted = _agg_df(df, sample_set, expression_col)
                print(df_pivoted.show())

                (df_pivoted.join(df_junction, sf.col('kmer') == sf.col('kmer_junction'), 'left_outer').
                 drop('kmer_junction').
                 coalesce(partitions).
                 write.mode('overwrite').
                 parquet(os.path.join(output_path))
                 )


def _agg_df(df, sample_set, agg_col):
    return (df.filter( f'max_{agg_col}' + '> 0' ).groupby('kmer').
            pivot('sample', values=list(sample_set)).
            agg(sf.max(f'max_{agg_col}').alias(f'max_{agg_col}'))
            )

def collapse_values(value):
    return np.nanmax(value.split('/'))



def mode_crosscohort(arg): #TODO one spark session or NOT?
    junction_field = 'junction_count'#'junction_expr'
    segment_field = 'expr'#'segment_expr'
    agg_fields_combine = [junction_field, segment_field, 'is_cross_junction']
    expression_fields = [junction_field, segment_field]
    grp_cols = 'kmer'
    if arg.remove_bg:
        backgr_suffix = "_no_back"
    else:
        backgr_suffix = ''
    if arg.compressed_inputs:
        compres_suffix = '.gz'

    if not arg.skip_filegrouping:
        logging.info(">>> Reorganize folder structure")
    paths_to_combine = organize_inputs(arg.input_dir, backgr_suffix, arg.kmer, arg.samples, arg.mutation_modes, compres_suffix, arg.skip_filegrouping)
    path_crosssamples = []
    logging.info(">>> Combine mutation modes per sample")
    for kmer_files, combined_path in paths_to_combine:
        run_combine(inf=kmer_files,
                    outf=combined_path,
                    partition_field=grp_cols,
                    grp_fields=grp_cols,
                    expression_fields=expression_fields,
                    max_aggregation_fields=agg_fields_combine,
                    cores = arg.cores,
                    memory_per_core=arg.mem_per_core,
                    partitions=1) # arg.cores * 10)
        path_crosssamples.append(combined_path)

    logging.info(">>> Combine kmer count info in kmers x sample matrix")
    run_pivot(input_=path_crosssamples,
              output_dir=arg.output_dir,
              output_suffix=arg.output_suffix, # os.path.join(arg.output_dir, "crosssamples_" + agg_fields_pivot + arg.output_suffix + '.pq'),
              sample_list=arg.samples,
              aggregation_col=expression_fields,
              partitions=1, #arg.cores * 10,
              cores=arg.cores,
              memory_per_core=arg.mem_per_core)




# All edge ans segment matrixes
#filter zeros
# add split '/'
# rewrite test
