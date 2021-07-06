###Line profiler

from immunopepper.utils import get_exon_expr
from line_profiler import LineProfiler
import pickle

out_dir = '/Users/laurieprelot/Documents/Projects/tmp_kmer/ERR2130621'
gene = pickle.load(open('/Users/laurieprelot/Documents/Projects/tmp_kmer/ERR2130621/sample_gene.pkl', 'rb' ))
vstart = 4493099
vstop = 4493406
countinfo = pickle.load(open('/Users/laurieprelot/Documents/Projects/tmp_kmer/ERR2130621/countinfo.pkl', 'rb' ))
Idx = pickle.load(open('/Users/laurieprelot/Documents/Projects/tmp_kmer/ERR2130621/idx.pkl', 'rb' ))
seg_counts = None



lp = LineProfiler()
lp_wrapper = lp(get_exon_expr)
lp_wrapper(gene, vstart, vstop, countinfo, Idx, seg_counts=None)
lp.print_stats()



#!kernprof -l -v target_function