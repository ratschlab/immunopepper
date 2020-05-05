import gzip
import logging
import pickle
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


def gz_and_normal_open(file_path,mode='r'):
    if file_path.endswith('.gz'):
        mode += 't'
        file_fp = gzip.open(file_path, mode)
    else:
        file_fp = open(file_path, mode)

    return file_fp


def write_list(fp, _list):
    fp.writelines([l + '\n' for l in _list])


def _convert_list_to_str(_list):
    remove_none_list = filter(lambda x:x is not None, _list)
    return ';'.join([str(_item) for _item in remove_none_list])

def convert_namedtuple_to_str(_namedtuple, field_list = None, sep = '\t'):

    if field_list is None:
        field_list = _namedtuple._fields
    line = ''
    for field in field_list:
        if field == 'new_line':
            line = line.strip() + '\n'
            continue
        if not hasattr(_namedtuple, field):
            logging.error('Namedtuple %s ' % str(_namedtuple) + ' does not have a field: ' + field)
            sys.exit(1)
        item = getattr(_namedtuple, field)
        if isinstance(item, (list, tuple)):
            line += _convert_list_to_str(item)+sep
        else:
            line += str(item)+sep
    return line[:-1] # remove the last '\t'


def write_namedtuple_list(fp, namedtuple_list, field_list):
    """ Write namedtuple_list to the given file pointer"""
    fp.writelines(convert_namedtuple_to_str(_namedtuple, field_list) + '\n' for _namedtuple in namedtuple_list)

def write_gene_expr(fp, gene_expr_tuple_list):
    header_line = 'gene\texpr\n'
    fp.write(header_line)
    gene_expr_str_list = [ gene_expr_tuple[0]+'\t'+str(gene_expr_tuple[1]) for gene_expr_tuple in gene_expr_tuple_list]
    write_list(fp,gene_expr_str_list)

def codeUTF8(s):
    return s.encode('utf-8')

def decodeUTF8(s):
    if not hasattr(s, 'decode'):
        return s
    return s.decode('utf-8')


