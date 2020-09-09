from collections import defaultdict
import glob
import logging
import os
import pathlib
from pyspark.sql import functions as sf
import shutil
import sys
from typing import Iterable

from .config import create_spark_session
from .config import create_spark_session_from_config
from .config import default_spark_config


def organize_inputs(input_dir, backgr_suffix, kmer, samples, mutation_modes, compres_suffix):
    for sample in samples:
        basedir = os.path.join(input_dir, sample)
        pathlib.Path(os.path.join(basedir, "junction_kmers")).mkdir(exist_ok= True, parents= True)
        jct_inputs = ['{}_junction_{}mer_{}.pq{}'.format(mutation, kmer, backgr_suffix, compres_suffix)
                      for mutation in mutation_modes]
        for input_file in jct_inputs:
            if not os.path.isfile(input_file):
                logging.error("Error: input file {} does not exist".format(input_file))
                sys.exit(1)
            shutil.move(input_file, os.path.join(basedir, "junction_kmers"))


def run_combine(inf: Iterable[str],
                outf: str,
                partition_field: str,
                grp_fields: Iterable[str],
                max_aggregation_fields: Iterable[str],
                cores: int, memory_per_core:int, partitions: int):
    spark_cfg = default_spark_config(cores, memory_per_core)

    with create_spark_session_from_config(spark_cfg) as spark:
        df = spark.read.parquet(*inf)
        agg_cols = [sf.max(c).alias('max_' + c) for c in max_aggregation_fields]

        df_g = (df.repartition(partitions, partition_field).
                        groupBy(grp_fields).
                        agg(*agg_cols)
               )

        df_g.write.mode('overwrite').parquet(outf)


def run_pivot(input_, output, sample_list, aggregation_col, partitions, cores, memory_per_core):
    #     part_id = int(os.path.basename(output).split('-')[1])

    #     print(input_)
    #     paths = input_

    #     print(f'checking {paths}')
    #     sample_set = set(sample_list.split(','))

    #     part_cnt = partitions
    #     print(f'input partitions {part_cnt}')
    sample_set = set(sample_list)
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

            print(aggregation_col)
            df_pivoted = _agg_df(df, sample_set, aggregation_col)
            print(df_pivoted.show())

            (df_pivoted.join(df_junction, sf.col('kmer') == sf.col('kmer_junction'), 'left_outer').
             drop('kmer_junction').
             coalesce(partitions).
             write.mode('overwrite').
             parquet(os.path.join(output))
             )


def _agg_df(df, sample_set, agg_col):
    return (df.groupby('kmer').
            pivot('sample', values=list(sample_set)).
            agg(sf.max(f'max_{agg_col}').alias(f'max_{agg_col}'))
            )


input_folder = ['/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/output/peptides_ccell_rerun_200707/15K_c5_allg_batch10_a18874e/TCGA-AR-A1AP-01A-11/test_pq/dummy']
output_file = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/output/peptides_ccell_rerun_200707/15K_c5_allg_batch10_a18874e/TCGA-AR-A1AP-01A-11/test_pq/dummy/combined2_pq.gz'

input_pairs = input_folder
agg_fields = ['junction_expr', 'segment_expr', 'is_cross_junction']
grp_cols = 'kmer'
run_combine(input_pairs, output_file, grp_cols, grp_cols, agg_fields, 1, 5000, 3)


agg_fields = 'junction_expr' # Or segment expression
#grp_cols = 'kmer'
sample_list = ['dummy1', 'dummy2']
partitions = 3
cores = 1
memory_per_core = 5000
input_ = ['/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/output/peptides_ccell_rerun_200707/15K_c5_allg_batch10_a18874e/TCGA-AR-A1AP-01A-11/test_pq/dummy/combined_pq.gz', # fakesample
          '/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/output/peptides_ccell_rerun_200707/15K_c5_allg_batch10_a18874e/TCGA-AR-A1AP-01A-11/test_pq/dummy2/combined_pq.gz']

output = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/output/peptides_ccell_rerun_200707/15K_c5_allg_batch10_a18874e/TCGA-AR-A1AP-01A-11/test_pq/pivot'
run_pivot(input_, output, sample_list, agg_fields, partitions, cores, memory_per_core)
