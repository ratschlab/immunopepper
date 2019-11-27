import psutil
import os
import logging

def print_memory_diags(disable_print=False):
    """
    Print memory diagnostics including the active resident set size
    """
    process = psutil.Process(os.getpid())
    memory = process.memory_info().rss/1000000000.0
    if not disable_print:
        logging.info('\tMemory usage: {:.3f} GB'.format(memory))
    return memory




