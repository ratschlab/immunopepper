# python library
import sys
from collections import namedtuple
# external library
import numpy as np
import h5py
import scipy as sp

# immuno module
from immuno_print import print_memory_diags
from utils import to_adj_succ_list,find_overlapping_cds_simple,attribute_list_to_dict,leq_strand,encode_chromosome

# Pre-process gene structures to aid fast generation of cross-junction peptides
# get gene.vertex_succ_list
# sort the vertex
# get gene.splicegraph.reading_frames

def genes_preprocess(genes, gene_cds_begin_dict):
    for gene_idx in range(genes.shape[0]):
        if gene_idx > 0 and gene_idx % 100 == 0:
            sys.stdout.write('.')
            if gene_idx % 1000 == 0:
                sys.stdout.write('%i/%i\n' % (gene_idx, genes.shape[0]))
            sys.stdout.flush()
        gene = genes[gene_idx]
        gene.from_sparse()
        assert (gene.strand in ["+", "-"])
        assert (len(gene.transcripts) == len(gene.exons))

        # Ignore genes that have no CDS annotated...
        if gene.name not in gene_cds_begin_dict.keys():
            continue

        gene.nvertices = gene.splicegraph.vertices.shape[1]
        gene.vertex_succ_list = to_adj_succ_list(gene.splicegraph.edges, gene.splicegraph.vertices, gene.strand)
        if gene.strand == "+":
            gene.vertex_order = np.argsort(gene.splicegraph.vertices[0, :])
        else:  # gene.strand=="-"
            gene.vertex_order = np.argsort(gene.splicegraph.vertices[1, :])[::-1]

        # get the reading_frames
        gene.splicegraph.reading_frames = {}
        gene.vertex_len_dict = {}
        for idx in gene.vertex_order:
            gene.splicegraph.reading_frames[idx] = set()
            v_start = gene.splicegraph.vertices[0, idx]
            v_stop = gene.splicegraph.vertices[1, idx]
            cds_begins = find_overlapping_cds_simple(v_start, v_stop, gene_cds_begin_dict[gene.name], gene.strand)
            gene.vertex_len_dict[idx] = v_stop - v_start
            if gene.vertex_len_dict[idx] < 3:
                continue
            # Initialize reading regions from the CDS transcript annotations
            for cds_begin in cds_begins:
                line_elems = cds_begin[2]
                cds_strand = line_elems[6]
                assert (cds_strand == gene.strand)
                cds_phase = int(line_elems[7])
                cds_left = int(line_elems[3])-1
                cds_right = int(line_elems[4])

                #TODO: need to remove the redundance of (cds_start, cds_stop, item)
                # It's ugly and not necessary
                if gene.strand == "-":
                    cds_right_modi = cds_right - cds_phase
                    cds_left_modi = v_start
                    n_trailing_bases = cds_right_modi - cds_left_modi
                else:
                    cds_left_modi = cds_left + cds_phase
                    cds_right_modi = v_stop
                    n_trailing_bases = cds_right_modi - cds_left_modi

                if n_trailing_bases < 3:
                    continue
                read_phase = n_trailing_bases % 3
                gene.splicegraph.reading_frames[idx].add((cds_left_modi, cds_right_modi, read_phase))

        gene.to_sparse()

# Pre-processed the annotation file and builds a lookup structure that can be used to retrieve
# the CDS start positions when looping over the genes
##TODO: do we really need so many dictionary?
def preprocess_ann(ann_path):

    Table = namedtuple('Table',['ts_to_cds','gene_to_ts'])
    transcript_to_gene_dict = {}    # transcript -> gene id
    gene_to_transcript_dict = {}    # gene_id -> list of transcripts
    transcript_to_cds_dict = {}     # transcript -> list of CDS exons
    transcript_cds_begin_dict = {}  # transcript -> first exon of the CDS
    gene_cds_begin_dict = {}        # gene -> list of first CDS exons

    # collect information from annotation file
    for line in open(ann_path, 'r'):
        item = line.split('\t')
        feature_type = item[2]
        attribute_list = item[-1].split('; ')
        attribute_dict = attribute_list_to_dict(attribute_list)

        # store relationship between gene ID and its transcript IDs
        if feature_type in ['transcript', 'mRNA']:
            gene_id = attribute_dict['gene_id']
            transcript_id = attribute_dict['transcript_id']
            assert (transcript_id not in transcript_to_gene_dict)
            transcript_to_gene_dict[transcript_id] = gene_id

            try:
                gene_to_transcript_dict[gene_id].add(transcript_id)
            except KeyError:
                gene_to_transcript_dict[gene_id] = set([transcript_id])

        # Todo python is 0-based while gene annotation file(.gtf, .vcf, .maf) is one based
        elif feature_type == "CDS":
            parent_ts = attribute_dict['transcript_id']
            strand_mode = item[6]
            cds_left = int(item[3])-1
            cds_right = int(item[4])-1
            frameshift = int(item[7])
            try:
                transcript_to_cds_dict[parent_ts].append((cds_left, cds_right, frameshift))
            except KeyError:
                transcript_to_cds_dict[parent_ts] = [(cds_left, cds_right, frameshift)]
            if strand_mode == "+" :
                cds_start, cds_stop = cds_left, cds_right
            else:
                cds_start, cds_stop = cds_right, cds_left

            # we only consider the start of the whole CoDing Segment
            if parent_ts not in transcript_cds_begin_dict or \
               leq_strand(cds_start, transcript_cds_begin_dict[parent_ts][0], strand_mode):
                transcript_cds_begin_dict[parent_ts] = (cds_start, cds_stop, item)

    # collect first CDS exons for all transcripts of a gene
    for ts_key in transcript_to_gene_dict:
        target_gene = transcript_to_gene_dict[ts_key]
        if target_gene not in gene_cds_begin_dict:
            gene_cds_begin_dict[target_gene] = []
        if ts_key in transcript_cds_begin_dict:
            gene_cds_begin_dict[target_gene].append(transcript_cds_begin_dict[ts_key])

    # sort list of CDS exons per transcript
    for ts_key in transcript_to_cds_dict:
        transcript_to_cds_dict[ts_key] = sorted(transcript_to_cds_dict[ts_key], key=lambda coordpair: coordpair[0])

    table = Table(transcript_to_cds_dict, gene_to_transcript_dict)
    return gene_cds_begin_dict, table


def search_edge_metadata_segmentgraph(gene, sorted_pos, edges, Idx):
    ''' Gives the ordered edge coordinates of the edge, return expression information of the edge'''
    gene_name = gene.name
    segmentgraph = gene.segmentgraph
    edge_idxs = edges.lookup_table[gene_name]

    a = sp.where(segmentgraph.segments[1, :] == sorted_pos[1])[0]
    b = sp.where(segmentgraph.segments[0, :] == sorted_pos[2])[0]
    if a < b:
        idx = sp.ravel_multi_index([a, b], segmentgraph.seg_edges.shape)
    else:
        idx = sp.ravel_multi_index([b, a], segmentgraph.seg_edges.shape)
    cidxs = list(filter(lambda elem: elem[1] == idx, edge_idxs))
    assert (len(cidxs) == 1)
    cidx = cidxs[0]
    count = edges.expr[cidx[0], Idx.sample]
    return count


## Given a segment graph and a 1-based genomic coordinate position of a somatic mutation, returns the expression
#  level of the one segment that contains the position. Performs a linear search in the segments contained in
#  the graph. Could be optimized to do a binary search, let's see if it is a bottleneck.
def search_metadata_segmentgraph(gene, gen_coord, seg_lookup_table, strain_idx_table, segment_expr_info, donor_id,
                                 sample_suffix=None):
    gene_name = gene.name
    segment_graph = gene.segmentgraph
    seg_idxs = seg_lookup_table[gene_name]
    col_idx = strain_idx_table[donor_id + sample_suffix]
    assert (len(seg_idxs) == segment_graph.segments.shape[1])

    # We use that in a segment graph, the vertices are not overlapping, so once we have found the segment,
    # we can look-up its expression value and return directly.
    for per_gene_idx, global_seg_idx in enumerate(seg_idxs):
        lbound = segment_graph.segments[0, per_gene_idx] + 1
        rbound = segment_graph.segments[1, per_gene_idx]

        if gen_coord >= lbound and gen_coord <= rbound:
            gene_expr_entry = float(segment_expr_info[global_seg_idx, col_idx])
            return gene_expr_entry, per_gene_idx

    assert (False)
    return np.nan

## Constructs look-up tables to lookup the correct column in matrix for a particular donor and
#  splits the global segments matrix into smaller chunks for each geneID to enable faster look-up
#  during peptide emission.
def parse_gene_metadata_info(h5f, donor_list):
    strain_expr_info = h5f["/strains"]
    segment_expr_info = h5f["/segments"]
    edge_expr_info = h5f["/edges"]
    edge_idx_info = h5f["/edge_idx"]
    assert (strain_expr_info.size == segment_expr_info.shape[1])
    strain_idx_table = {}

    Segments = namedtuple('Segments',['expr','lookup_table'])
    Edges = namedtuple('Edges',['expr','lookup_table'])

    #TODO: make it clear how strain_id come from in h5f file
    for strain_idx in np.arange(strain_expr_info.size):
        strain_id = strain_expr_info[strain_idx]
        #strain_id = '-'.join(strain_id.split('.')[0].split('-')[:3])
        strain_id = 'test1'
        if strain_id in donor_list:
            strain_idx_table[strain_id] = strain_idx

    gene_names = h5f["/gene_names"]
    gene_ids_segs = h5f["/gene_ids_segs"]
    gene_ids_edges = h5f["/gene_ids_edges"]
    assert (gene_ids_segs.size == segment_expr_info.shape[0])
    assert (gene_ids_edges.size == edge_expr_info.shape[0])
    seg_lookup_table = {}
    edge_lookup_table = {}

    #print(gene_ids_segs.shape)
    for seg_idx in np.arange(gene_ids_segs.shape[0]):
        gene_id = gene_names[gene_ids_segs[seg_idx, 0]][0]  # no [0] before
        if gene_id not in seg_lookup_table:
            seg_lookup_table[gene_id] = []
        seg_lookup_table[gene_id].append(seg_idx)

    for edge_idx in np.arange(gene_ids_edges.shape[0]):
        gene_id = gene_names[gene_ids_edges[edge_idx, 0]][0]  # no [0] before
        if gene_id not in edge_lookup_table:
            edge_lookup_table[gene_id] = []
        edge_lookup_table[gene_id].append((edge_idx, edge_idx_info[edge_idx]))

    segments = Segments(segment_expr_info,seg_lookup_table)
    edges = Edges(edge_expr_info,edge_lookup_table)

    return segments, edges, strain_idx_table


def parse_mutation_from_vcf(vcf_path):
    """
    Extract germline mutation information from given vcf file.

    Parameters
    ----------
    vcf_path: str, vcf file path

    Returns
    -------
    mut_dict: with key (sample, chromo) and values (var_dict)

    """
    print(vcf_path)
    f = open(vcf_path,'r')
    lines = f.readlines()
    mutation_dic = {}
    for line in lines:
        if line.strip()[:2] == '##':  # annotation line
            continue
        if line.strip()[0] == '#':  # head line
            fields = line.strip().split('\t')
            sample_list = fields[9:]
            continue
        items = line.strip().split('\t')
        var_dict = {}
        chr = items[0]
        pos= int(items[1])-1
        var_dict['ref_base'] = items[3]
        var_dict['mut_base'] = items[4]
        var_dict['qual'] = items[5]
        var_dict['filter'] = items[6]
        info = items[7]
        info_items =info.split(';')
        for info_item in info_items:
            key, value = info_item.split('=')
            if key == 'VT':
                var_dict['varianttype'] = value
                break
        if var_dict['varianttype'] == 'SNP':  # only consider snp for now
            for sample_id in sample_list:
                if (sample_id,chr) in mutation_dic.keys():
                    mutation_dic[(sample_id,chr)][int(pos)] = var_dict
                else:
                    mutation_dic[(sample_id,chr)] = {}
                    mutation_dic[(sample_id, chr)][int(pos)] = var_dict
    return mutation_dic

def parse_mutation_from_vcf_h5(h5_vcf_path,sample_list):
    """
    Extract germline mutation information from given vcf h5py file.

    Parameters
    ----------
    h5_vcf_path: str, vcf file path
    sample_list: list of str, list for sample name

    Returns
    -------
    mut_dict: with key (sample, chromo) and values (var_dict)

    """
    a = h5py.File(h5_vcf_path,'r')
    mut_dict = {}
    for sample in sample_list:
        col_id = [i for (i, item) in enumerate(a['gtid']) if item.startswith(sample)][0]
        row_id = sp.where(a['gt'][:,col_id]==0)[0]
        for irow in row_id:
            chromo = encode_chromosome(a['pos'][irow,0])
            pos = a['pos'][irow,1]-1
            mut_base = a['allele_alt'][irow]
            ref_base = a['allele_ref'][irow]
            var_dict = {"mut_base":mut_base,"ref_base":ref_base}
            if (sample,chromo)  in mut_dict:
                mut_dict[(sample,chromo)][pos] = var_dict
            else:
                mut_dict[(sample,chromo)] = {}
                mut_dict[(sample,chromo)][pos] = var_dict
    return mut_dict

def parse_mutation_from_maf(maf_path):
    '''
    Extract somatic mutation information from given maf file.

    Parameters
    ----------
    maf_path: str, maf file path

    Returns
    -------
    mut_dict: with key (sample, chromo) and values (var_dict)

    '''
    f = open(maf_path)
    lines = f.readlines()
    mutation_dic = {}
    for i,line in enumerate(lines[1:]):
        print(i)
        items = line.strip().split('\t')
        if items[9] == 'SNP':  # only consider snp
            sample_id = '-'.join(items[15].split('-')[:3])
            chr = items[4]
            pos = int(items[5])-1
            var_dict = {}
            var_dict['ref_base'] = items[10]
            var_dict['mut_base'] = items[12]
            var_dict['strand'] = items[7]
            var_dict['variant_Classification'] = items[8]
            var_dict['variant_Type'] = items[9]
            if (sample_id,chr) in mutation_dic.keys():
                mutation_dic[((sample_id,chr))][int(pos)] = var_dict
            else:
                mutation_dic[((sample_id, chr))] = {}
                mutation_dic[((sample_id,chr))][int(pos)] = var_dict
    return mutation_dic

def parse_junction_meta_info(h5f_path):
    """
    Extract introns of interest from given h5py file
    Parameters
    ----------
    h5f_path: str, h5py file path

    Returns
    -------
    junction_dict: dict, key (chromosome id), value (set of coordinate pairs)

    """
    if h5f_path is None:
        return None
    else:
        h5f = h5py.File(h5f_path)
        chrms = h5f['chrms'][:]
        pos = h5f['pos'][:].astype('str')
        strand = h5f['strand'][:]
        junction_dict = {}

        for i,ichr in enumerate(chrms):
            try:
                junction_dict[ichr].add(':'.join([pos[i, 0], pos[i, 1], strand[i]]))
            except KeyError:
                junction_dict[ichr] = set([':'.join([pos[i, 0], pos[i, 1], strand[i]])])
    return junction_dict









