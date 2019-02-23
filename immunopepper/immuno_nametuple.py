from collections import namedtuple

Output_peptide = namedtuple('Output_peptide', ['output_id','vertex_pair_id','peptide'])
Output_metadata = namedtuple('Output_metadata', ['output_id', 'read_frame', 'gene_name',	'gene_chr',
                                                 'gene_strand',	'mutation_mode','peptide_weight', 'peptide_annotated'
                                                 'junction_annotated',	'has_stop_codon',
                                                 'is_in_junction_list',	'is_isolated',
                                                 'variant_comb',	'variant_seg_expr',
                                                 'exons_coor',	'vertex_idx',	'junction_expr',
                                                 'segment_expr'])
Output_background = namedtuple('Output_background', ['trans_id', 'peptide'])
#expr = namedtuple('expr', [])



