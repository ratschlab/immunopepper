from pyspark.sql import SparkSession
from pyspark.conf import SparkConf
import tblib.pickling_support

tblib.pickling_support.install()


def default_spark_config(cores: int, memory_per_executor_mb: int, parallelism: int,
                         tmp_dir: str = '', extra_java_options: str = '', enable_arrow: bool = True,
                         use_utc: bool = False) -> SparkConf:
    """
    Creates a default configuration for running Apache Spark.
    See also https://spark.apache.org/docs/latest/configuration.html for more information.
    :param cores: number of executors (workers)
    :param memory_per_executor_mb: memory per executor [MB]
    :param tmp_dir: "scratch" space for spark, typically /tmp by default
    :param extra_java_options: extra parameters for the JVM
    :param enable_arrow: see https://spark.apache.org/docs/2.3.0/sql-programming-guide.html#pyspark-usage-guide-for-pandas-with-apache-arrow ,
           requires pyarrow to be installed
    :return: SparkConf instance
    """
    driver_mem = int(0.75 * cores * memory_per_executor_mb)
    memory_per_executor_mb = int(memory_per_executor_mb * 0.8)
    print("driver_mem", driver_mem)
    print("memory_per_executor_mb 80%", memory_per_executor_mb)
    print("parallelism_", parallelism)
    print("shuffle_partitions", parallelism)
    print("permsize", "1024M")

    cfg = SparkConf()

    if tmp_dir:
        cfg.set("spark.local.dir", tmp_dir)

    # TODO set as parameter
    java_options = str(extra_java_options)
    #java_options = java_options + "-verbose:gc -XX:+PrintGCDetails -XX:+PrintGCTimeStamps -XX:+PrintFlagsFinal"
    java_options = java_options + " -XX:ThreadStackSize=81920"

    if use_utc:
        # avoiding trouble with JDBC and timestamps
        java_options = "-Duser.timezone=UTC " + java_options

    return (cfg.set("spark.driver.memory", "{}m".format(driver_mem)).
            set("spark.executor.memory", "{}m".format(memory_per_executor_mb)).
            set("spark.driver.extraJavaOptions", java_options).
            set("spark.master", "local[{}]".format(cores)).
            # set("spark.executor.cores", int(np.floor(cores / core_per_exec))).
            # set("spark.jars", jar_paths).
            set("spark.sql.execution.arrow.pyspark.enabled", str(enable_arrow)).  # TODO set as parameter
            set("spark.sql.debug.maxToStringFields", 11000).
            set("spark.executor.heartbeatInterval", 10000).
            set("spark.network.timeout", 1000000).
            set("spark.serializer", "org.apache.spark.serializer.KryoSerializer").
            # set("spark.driver.bindAddress", "192.168.0.13").
            set("spark.default.parallelism", parallelism).
            set("spark.sql.shuffle.partitions", parallelism).
            set("spark.driver.maxResultSize", "0")  # unlimited
            # TODO remove the personal IP address
            )


def create_spark_session(cores: int, memory_per_executor_mb: int) -> SparkSession:
    '''
    Creates a local spark session with a given number of executors and memory
    :param cores: number of executors (workers)
    :param memory_per_executor_mb: memory per executor [MB]. A default overhead for the driver is added
    :return: SparkSession instance
    '''
    return create_spark_session_from_config(default_spark_config(cores, memory_per_executor_mb))


def create_spark_session_from_config(cfg: SparkConf) -> SparkSession:
    return (SparkSession.
            builder.
            config(conf=cfg).
            getOrCreate())
