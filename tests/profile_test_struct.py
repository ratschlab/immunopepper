import os
import pstats

tmpdir = '/Users/laurieprelot/Documents/Projects/tmp_kmer'
pr = pstats.Stats(os.path.join(tmpdir, 'cProfile.pstats'))
pr.sort_stats('tottime').print_stats('immunopepper')
