import psutil
import os

## Overview of somatic mutation information
def print_som_info(som_table):
    for k, v in som_table.iteritems():
        print("Somatic mutations for donor: {}".format(k))

        print("Number of mutations: {}".format(len(v)))

def print_memory_diags():
    ''' Print memory diagnostics including the active resident set size'''
    process = psutil.Process(os.getpid())
    print('\tMemory usage: {:.3f} GB\n'.format(process.memory_info().rss/1000000000.0))




