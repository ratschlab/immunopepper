from collections import namedtuple


"""
Coord namedtuple
start and stop position of the junction pairs (two exons). If it only consists of one exon,
start_v2 and stop_v2 is NOT_EXIST.
"""
try:
    Coord = namedtuple('Coord', ['start_v1', 'stop_v1', 'start_v2', 'stop_v2', 'start_v3', 'stop_v3'], defaults=(None,) * 2)
except:
    Coord = namedtuple('Coord', ['start_v1', 'stop_v1', 'start_v2', 'stop_v2', 'start_v3', 'stop_v3'])
    Coord.__new__.__defaults__ = (None,) * 2

"""
Output_junc_peptide namedtuple
- output_id: (gene_id).(junction_id). gene_id is the index of given gene in the splicegraph array.
    junction_id is the index of given junction pair in all junction pair (in descending or ascending order)
- id: (gene_name)_(first_vertex_id)_(second_vertex_id). Detail information of output_id. 
    We can know clearly from which gene and which vertex pair the given peptide is translated.
- peptide: (peptide_string). The peptide translated from junction pairs.
- exons_coor: Coord namedtuple
"""
OutputJuncPeptide = namedtuple('OutputJuncPeptide', ['output_id','id','peptide','exons_coor','junction_count'])


"""
Output_metadata namedtuple. 
- output_id: the same with that in Output_junc_peptide
- read_frame: int (0,1,2). The number of base left to the next junction pair. 
- gene_name: str. The name of Gene.
- gene_chr: str. The Chromosome id where the gene is located.
- gene_strand: str ('+', '_'). The strand of gene.
- mutation_mode: str ('ref', 'somatic', 'germline', 'somatic_and_germline'). Mutation mode
- peptide_annotated: Boolean. Indicate if the junction peptide also appears in the background peptide.
- junction_peptided: Boolean. Indicate if the junction also appear in the input annotation file.
- has_stop_codon. Boolean. Indicate if there is stop codon in the junction pair.
- is_in_junction_list: Boolean. Indicate if the junction pair appear in the given junction whitelist
- is_isolated: Boolean. Indicate if the output peptide is actually translated from a single exon instead of two.
- variant_comb: shows the somatic variantion combination used in this line of output. seperate by ';'
    eg. 5;25 means the somatic mutation of position 5 and 25 take effect in this output.
- variant_seg_expr: shows the corresponding expression of segments where the corresponding somatic mutation is in.
    eg. 257.0;123.2 means the segment where the somatic mutation in position 5 is in has counts 257.0
- modified_exons_coor. Coord namedtuple. Shows exon coordination. Usually we have 4 number start_v1;stop_v1;start_v2;stop_v2. They 
    have already absorb reading frame so you can use the coord directly to generate the same output peptide.
- original_exons_coord. Coord namedtuple. Shows the original exon coordination.
- vertex_idx. shows the vertex id of the given junction. eg 5,6 means this junction pars consists of the fifth and
    sixth vertex.
- junction_expr. float. The expression of the junction.
- segment_expr. float. The weighted sum of segment expression. We split the junction into segments and compute the segment 
    expression with the length-weighted-sum expression.
"""
OutputMetadata = namedtuple('OutputMetadata', ['output_id', 'read_frame', 'gene_name', 'gene_chr',
                                                 'gene_strand',	'mutation_mode', 'peptide_annotated',
                                                 'junction_annotated',	'has_stop_codon',
                                                 'is_in_junction_list',	'is_isolated',
                                                 'variant_comb',	'variant_seg_expr',
                                                 'modified_exons_coord','original_exons_coord',	'vertex_idx',	'junction_expr',
                                                 'segment_expr'])

VertexPair = namedtuple('VertexPair', ['output_id', 'read_frame','has_stop_codon','modified_exons_coord','original_exons_coord','vertex_idxs','peptide_weight'])


"""
Output_backgrouond namedtuple.
- id: transcript name
- peptide: background peptide
"""
OutputBackground = namedtuple('OutputBackground', ['id', 'peptide'])


"""
Output_kmer namedtuple.
- kmer: output kmer
- id: transcript id(generated from background peptide) or gene_vertex_id (generated from concat peptide)
- expr: float. length-weighted sum of expression of the kmer
- is_cross_junction: boolen. indicate if the kmer spans over the cross junction
"""
OutputKmer= namedtuple('OutputKmer', ['kmer','id','expr','is_cross_junction','junction_count'])


"""
Peptide namedtuple
- mut: peptide translated from mutated dna string
- ref: peptide translated from original dna string
"""
Peptide = namedtuple('Peptide', ['mut', 'ref'])


"""
Flag namedtuple
- has_stop: Boolean variable that indicating if there is stop codon in the junction pair.
- is_isolated: Boolean variable indicatimng if the output peptide is actually translated from a single exon instead of two.
"""
Flag = namedtuple('Flag', ['has_stop', 'is_isolated'])


"""
Idx namedtuple
- gene: gene index in input splicegraph
- sample: sample index in count hdf5 mat.
"""
Idx = namedtuple('Idx', ['gene', 'sample'])


"""
Reading_frame_tuple namedtuple
- cds_left_modi: modified left cds coordinate. (modifies means read frame shift has already been considered)
- cds_right_modi: modified right cds coordinate
- read_phase: (0,1,2). the number of bases left for the next cds
"""
ReadingFrameTuple = namedtuple('ReadingFrameTuple',['cds_left_modi','cds_right_modi','read_phase'])


"""
GeneTable namedtuple
- gene_to_cds_begin: dict. (gene_name -> cds_begin)
- ts_to_cds: dict. (transcript_name -> List[Reading_frame_tuple])
- gene_to_ts: dict. (gene_name -> [transcript_name])
"""
GeneTable = namedtuple('GeneTable', ['gene_to_cds_begin', 'ts_to_cds', 'gene_to_ts'])


"""
CountInfo namedtuple
- sample_idx_dict: Dictionary mapping a sample name to an index
- gene_idx_dict: Dictionary mapping a gene name to an index
- gene_id_to_segrange: map gene ids to the appropriate segment range in count HDF5
- gene_id_to_edgerange: map gene ids to the appropriate edge range in count HDF5
- h5f: file name of the count hdf5 file
"""
CountInfo = namedtuple('CountInfo', ['sample_idx_dict', 'gene_idx_dict', 'gene_id_to_segrange', 'gene_id_to_edgerange', 'h5fname'])


"""
filepointer namedtuple
namedtuple that contain all the filepointers
- junction_peptide_fp:
- junction_meta_fp:
- background_peptide_fp:
"""
Filepointer = namedtuple('Filepointer',['junction_peptide_fp','junction_meta_fp','background_peptide_fp'])


"""
Mutation namedtuple
namedtuple that contains all mutation information
- mode: 
- germline_mutation_dict:
- somatic_mutation_dict:
"""
Mutation = namedtuple('Mutation', ['mode', 'germline_mutation_dict', 'somatic_mutation_dict'])

GeneInfo = namedtuple('GeneInfo', ['vertex_succ_list', 'vertex_order', 'reading_frames', 'vertex_len_dict', 'nvertices'])
