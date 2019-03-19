import psutil
import os


def print_memory_diags():
    """
    Print memory diagnostics including the active resident set size
    """
    process = psutil.Process(os.getpid())
    print('\tMemory usage: {:.3f} GB\n'.format(process.memory_info().rss/1000000000.0))




