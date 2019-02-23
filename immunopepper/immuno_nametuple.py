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
Peptide = namedtuple('Peptide', ['mut', 'ref'])
Coord = namedtuple('Coord', ['start_v1', 'stop_v1', 'start_v2', 'stop_v2'])
Flag = namedtuple('Flag', ['has_stop', 'is_isolated'])
Idx = namedtuple('Idx', ['gene', 'sample'])

#expr = namedtuple('expr', [])



