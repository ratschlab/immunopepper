import multiprocessing as mp
from pyspark.sql import SparkSession
from pyspark.conf import SparkConf
import sys
import threading
import tblib.pickling_support


tblib.pickling_support.install()




class ExceptionWrapper(object):

    def __init__(self, ee):
        self.ee = ee
        __, __, self.tb = sys.exc_info()

    def re_raise(self):
        raise self.ee.with_traceback(self.tb)



class MaxQueuePool:
    def __init__(self, *args, **kwargs):
        self.pool = mp.Pool( *args, **kwargs)

    def submit(self, function, my_args):
        """Submits a new task to the pool, blocks if Pool queue is full."""
        res = self.pool.apply_async(function, args=my_args) #, callback=self.task_done)
        error_message = res.get()
        if isinstance(error_message, ExceptionWrapper):
            error_message.re_raise()

    def terminate(self):
        self.pool.close()
        self.pool.join()


def default_spark_config(cores: int, memory_per_executor: int, driver_overhead: int = 2000,
                         tmp_dir: str = '', extra_java_options: str = '', enable_arrow: bool = True, use_utc: bool= False) -> SparkConf:
    '''
    See also https://spark.apache.org/docs/latest/configuration.html for more information
    about the semantics of the configurations
    :param cores: number of executors (workers)
    :param memory_per_executor: memory per executor [MB]
    :param driver_overhead: Memory to allocate for the driver [MB], excluding exector memory
    :param tmp_dir: "scratch" space for spark, typically /tmp by default
    :param extra_java_options: extra parameters for the JVM
    :param enable_arrow: see https://spark.apache.org/docs/2.3.0/sql-programming-guide.html#pyspark-usage-guide-for-pandas-with-apache-arrow , requires pyarrow to be installed
    :return: SparkConf instance
    '''
    driver_mem = cores * memory_per_executor + driver_overhead

    cfg = SparkConf()

    if tmp_dir:
        cfg.set("spark.local.dir", tmp_dir)

    java_options = str(extra_java_options)

    if use_utc:
        # avoiding trouble with JDBC and timestamps
        java_options = "-Duser.timezone=UTC " + java_options

    return (cfg.set("spark.driver.memory", "{}m".format(driver_mem)).
            set("spark.executor.memory", "{}m".format(memory_per_executor)).
            set("spark.driver.extraJavaOptions", java_options).
            set("spark.master", "local[{}]".format(cores)).
           # set("spark.jars", jar_paths).
            set("spark.sql.execution.arrow.pyspark.enabled", str(enable_arrow))#.
            #set("spark.driver.bindAddress", "192.168.0.14") #TODO remove the personal IP address
            )


def create_spark_session(cores: int, memory_per_executor: int) -> SparkSession:
    '''
    Creates a local spark session with a given number of executors and memory
    :param cores: number of executors (workers)
    :param memory_per_executor: memory per executor [MB]. A default overhead for the driver is added
    :return: SparkSession instance
    '''
    return create_spark_session_from_config(default_spark_config(cores, memory_per_executor))


def create_spark_session_from_config(cfg: SparkConf) -> SparkSession:
    return (SparkSession.
             builder.
             config(conf=cfg).
             getOrCreate())
