import pickle
import sys

import spladder.classes.gene as gene
import spladder.classes.splicegraph as splicegraph
import spladder.classes.segmentgraph as segmentgraph
sys.modules['modules.classes.gene'] = gene
sys.modules['modules.classes.splicegraph'] = splicegraph
sys.modules['modules.classes.segmentgraph'] = segmentgraph

# TODO: this is a hack to get around module name changes spladder affecting pickle
def load_pickled_graph(f):

    return pickle.load(f)

