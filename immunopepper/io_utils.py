import pickle


# TODO: this is a hack to get around module name changes spladder affecting pickle
def load_pickled_graph(f):
    # following https://wiki.python.org/moin/UsingPickle/RenamingModules
    renametable = {
        'modules.classes.gene': 'spladder.classes.gene',
        'modules.classes.splicegraph': 'spladder.classes.splicegraph',
        'modules.classes.segmentgraph': 'spladder.classes.segmentgraph'
    }

    def mapname(name):
        if name in renametable:
            return renametable[name]
        return name

    def mapped_load_global(self):
        module = mapname(self.readline()[:-1])
        name = mapname(self.readline()[:-1])
        klass = self.find_class(module, name)
        self.append(klass)

    unpickler = pickle.Unpickler(f)

    unpickler.dispatch[pickle.GLOBAL] = mapped_load_global
    return unpickler.load()
