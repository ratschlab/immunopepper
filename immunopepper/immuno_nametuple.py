from collections import namedtuple

Output_junc_peptide = namedtuple('Output_junc_peptide', ['output_id','id','peptide'])
Output_metadata = namedtuple('Output_metadata', ['output_id', 'read_frame', 'gene_name',	'gene_chr',
                                                 'gene_strand',	'mutation_mode','peptide_weight', 'peptide_annotated',
                                                 'junction_annotated',	'has_stop_codon',
                                                 'is_in_junction_list',	'is_isolated',
                                                 'variant_comb',	'variant_seg_expr',
                                                 'exons_coor',	'vertex_idx',	'junction_expr',
                                                 'segment_expr'])
Output_background = namedtuple('Output_background', ['id', 'peptide'])
Output_kmer= namedtuple('Output_kmer', ['kmer','id','expr'])
Peptide = namedtuple('Peptide', ['mut', 'ref'])
Coord = namedtuple('Coord', ['start_v1', 'stop_v1', 'start_v2', 'stop_v2'])
Flag = namedtuple('Flag', ['has_stop', 'is_isolated'])
Idx = namedtuple('Idx', ['gene', 'sample'])
Reading_frame_tuple = namedtuple('reading_frame_tuple',['cds_left_modi','cds_right_modi','read_phase'])
GeneTable = namedtuple('GeneTable', ['gene_to_cds_begin', 'ts_to_cds', 'gene_to_ts'])
Segments = namedtuple('Segments', ['expr', 'lookup_table'])
Edges = namedtuple('Edges', ['expr', 'lookup_table'])
CountInfo = namedtuple('CountInfo', ['segments', 'edges', 'strain_idx_table'])




