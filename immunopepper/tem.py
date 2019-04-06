# load maf dict

from immunopepper.immuno_preprocess import  parse_gene_metadata_info
from immunopepper.immuno_mutation import get_sub_mutation_tuple
from immunopepper.utils import get_idx

import h5py
import cPickle
pickle_path = '/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018_alt_splice/spladder/genes_graph_conf2.merge_graphs.validated.pickle'
with open(pickle_path, 'rb') as graph_fp:
    (graph_data, graph_meta) = cPickle.load(graph_fp)  # both graph data and meta data
genes_preprocess(graph_data, genetable.gene_to_cds_begin)

# load genetable
from immunopepper.immuno_preprocess import preprocess_ann
ann_path = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/annotation/gencode.v19.annotation.hs37d5_chr.gff'
genetable = preprocess_ann(ann_path)

gene_id_dict = {}
for i,gene in enumerate(graph_data):
    gene_id_dict[gene.name] = i

# load sequence file
import Bio.SeqIO as BioIO
ref_path = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/sequence/genome.fa'
seq_dict = {}
interesting_chr = map(str, range(1, 23)) + ["X", "Y", "MT"]
print('Parsing genome sequence ...')
for record in BioIO.parse(ref_path, "fasta"):
    if record.id in interesting_chr:
        seq_dict[record.id] = str(record.seq).strip()


# load gene_graph
samples = ['TCGA-A2-A0D0']
from immunopepper.immuno_preprocess import parse_junction_meta_info
h5f = h5py.File('/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018_alt_splice/spladder/genes_graph_conf2.merge_graphs.validated.count.hdf5','r')
countinfo = parse_gene_metadata_info(h5f, samples)
edges, segments, sample_idx_table = countinfo.edges,countinfo.segments,countinfo.sample_idx_table


f = open('/cluster/work/grlab/projects/TCGA/immunopepper_rerun/maf.pickle')
maf_dict = cPickle.load(f)
from collections import namedtuple
Mutation = namedtuple('Mutation', ['mode', 'maf_dict', 'vcf_dict'])
mutation = Mutation('somatic',maf_dict,{})

gene_idx = 1221
gene=graph_data[gene_idx]
sub_mutation = get_sub_mutation_tuple(mutation,samples[0], gene.chr)
idx = get_idx(sample_idx_table, samples[0], gene_idx)


from immunopepper.utils import get_sub_mut_dna
from immunopepper.immuno_nametuple import Reading_frame_tuple
from immunopepper.utils import translate_dna_to_peptide
a = calculate_output_peptide(gene,seq_dict[gene.chr],idx,segments,edges,sub_mutation,genetable,debug=True)


from immunopepper.immuno_preprocess import preprocess_ann,genes_preprocess
import cPickle
ann_path = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/annotation/gencode.v19.annotation.hs37d5_chr.gff'
splice_path = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/graph/genes_graph_gtex_conf2.merge_graphs.pickle'
print('Building lookup structure ...')
genetable = preprocess_ann(ann_path)

# load splicegraph
print('Loading splice graph ...')
with open(splice_path, 'rb') as graph_fp:
    (graph_data, graph_meta) = cPickle.load(graph_fp)  # both graph data and meta data

genes_preprocess(graph_data, genetable.gene_to_cds_begin)

# check germ memory explosion
pos_start = np.min(gene.splicegraph.vertices)
pos_end = np.max(gene.splicegraph.vertices)
mut_dict =sub_dict
variant_pos_candi = [ipos for ipos in list(mut_dict.keys()) if ipos >= pos_start and ipos < pos_end]
print(len(variant_pos_candi))

num_v = [gene.splicegraph.vertices.shape[1] for gene in graph_data]
import numpy as np
large_graph_id = np.where(np.array(num_v)>5000)[0]
large_graph = graph_data[large_graph_id]

f = open('genes_graph_gtex_conf2.merge_graphs_large_graph.pickle','wb')
cPickle.dump((large_graph,graph_data))
