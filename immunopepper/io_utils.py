import gzip
import logging
import os
import pickle
import psutil
import sys

### intermediate fix to load pickle files stored under previous version
import spladder.classes.gene as gene
import spladder.classes.splicegraph as splicegraph
import spladder.classes.segmentgraph as segmentgraph
sys.modules['modules.classes.gene'] = gene
sys.modules['modules.classes.splicegraph'] = splicegraph
sys.modules['modules.classes.segmentgraph'] = segmentgraph
### end fix


def load_pickled_graph(f):
    return pickle.load(f)


def gz_and_normal_open(file_path):
    if file_path.endswith('.gz'):
        file_fp = gzip.open(file_path, 'wt')
    else:
        file_fp = open(file_path,'w')
    return file_fp


def print_memory_diags(disable_print=False):
    """
    Print memory diagnostics including the active resident set size
    """
    process = psutil.Process(os.getpid())
    memory = process.memory_info().rss/1000000000.0
    if not disable_print:
        logging.info('\tMemory usage: {:.3f} GB'.format(memory))
    return memory


