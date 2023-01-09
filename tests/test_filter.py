"""
Unit tests for imunopepper/filter.py
"""
import numpy as np
import pytest
from spladder.classes.gene import Gene
from spladder.classes.splicegraph import Splicegraph

import immunopepper.filter as filter
from immunopepper.namedtuples import Coord, OutputKmer, ReadingFrameTuple, VertexPair


class TestFilterRedundantJunctions:
    """ Tests for :meth:`filter.filter_redundant_junctions` """
    def test_no_redundant(self):
        read_frame = ReadingFrameTuple(cds_left_modi=4495135, cds_right_modi=4495155, read_phase=2, annotated_RF=True)
        # the second VertexPair is subsuming the first, because they both refer to the same read_phase (2), and the
        # same intron (4493490, 4495135) but the exon of the 2nd [4491713-4495155]
        # includes the coordinates of the first [4493099-4495155]
        vertex_pairs = [VertexPair(output_id='0:0',
                                   read_frame=read_frame, has_stop_codon=True,
                                   modified_exons_coord=Coord(start_v1=1, stop_v1=2, start_v2=3, stop_v2=4,
                                                              start_v3=None, stop_v3=None),
                                   original_exons_coord=None, vertex_idxs=[32, 14], peptide_weight='1.000'),
                        VertexPair(output_id='0:71',
                                   read_frame=read_frame, has_stop_codon=True,
                                   modified_exons_coord=Coord(start_v1=5, stop_v1=6, start_v2=7, stop_v2=8,
                                                              start_v3=None, stop_v3=None),
                                   original_exons_coord=None, vertex_idxs=[26, 5], peptide_weight='1.000')
                        ]
        assert vertex_pairs == filter.filter_redundant_junctions(vertex_pairs, '-')

    def test_one_redundant(self):
        read_frame = ReadingFrameTuple(cds_left_modi=4495135, cds_right_modi=4495155, read_phase=2, annotated_RF=True)
        # the second VertexPair is subsuming the first, because they both refer to the same read_phase (2), and the
        # same intron (4493490, 4495135) but the exon of the 2nd [4491713-4495155]
        # includes the coordinates of the first [4493099-4495155]
        vertex_pairs = [VertexPair(output_id='0:0',
                                   read_frame=read_frame, has_stop_codon=True,
                                   modified_exons_coord=Coord(start_v1=4495135, stop_v1=4495155, start_v2=4493099,
                                                              stop_v2=4493490, start_v3=None, stop_v3=None),
                                   original_exons_coord=None, vertex_idxs=[32, 14], peptide_weight='1.000'),
                        VertexPair(output_id='0:71',
                                   read_frame=read_frame, has_stop_codon=True,
                                   modified_exons_coord=Coord(start_v1=4495135, stop_v1=4495155, start_v2=4491713,
                                                              stop_v2=4493490, start_v3=None, stop_v3=None),
                                   original_exons_coord=None, vertex_idxs=[26, 5], peptide_weight='1.000')
                        ]
        filtered_vertex_pairs = filter.filter_redundant_junctions(vertex_pairs, '-')
        assert filtered_vertex_pairs == [vertex_pairs[1]]

    def test_identical_introns(self):
        read_frame = ReadingFrameTuple(cds_left_modi=4495135, cds_right_modi=4495155, read_phase=2, annotated_RF=True)
        # the second VertexPair is subsuming the first, because they both refer to the same read_phase (2), and the
        # same intron (4493490, 4495135) but the exon of the 2nd [4491713-4495155]
        # includes the coordinates of the first [4493099-4495155]
        vertex_pairs = [VertexPair(output_id='0:0',
                                   read_frame=read_frame, has_stop_codon=True,
                                   modified_exons_coord=Coord(start_v1=4495135, stop_v1=4495155, start_v2=4493099,
                                                              stop_v2=4493490, start_v3=None, stop_v3=None),
                                   original_exons_coord=None, vertex_idxs=[32, 14], peptide_weight='1.000'),
                        VertexPair(output_id='0:71',
                                   read_frame=read_frame, has_stop_codon=True,
                                   modified_exons_coord=Coord(start_v1=4495135, stop_v1=4495155, start_v2=4493099,
                                                              stop_v2=4493490, start_v3=None, stop_v3=None),
                                   original_exons_coord=None, vertex_idxs=[26, 5], peptide_weight='1.000')
                        ]
        filtered_vertex_pairs = filter.filter_redundant_junctions(vertex_pairs, '-')
        assert filtered_vertex_pairs == [vertex_pairs[1]]


class TestJunctionIsAnnotated:
    """ Tests for :meth:`filter.junction_is_annotated` """

    @pytest.fixture
    def gene(self):
        gene = Gene('ENSMUSG00000025902.13', 1000, 1200, 'chr1', '+')
        return gene

    @pytest.fixture
    def gene_to_transcript(self):
        return {
            'ENSMUSG00000025902.13': ['ENSMUST00000027035.9', 'ENSMUST00000195555.1', 'ENSMUST00000192650.5'],
            'ENSMUSG00000025903.14': ['ENSMUST00000027036.10', 'ENSMUST00000150971.7', 'ENSMUST00000119612.8']}

    @pytest.fixture
    def transcript_to_cds(self):
        return {'ENSMUST00000027035.9': [(4491718, 4492668, 2), (4493099, 4493406, 0)],
                'ENSMUST00000195555.1': [(4491718, 4492591, 0)],
                'ENSMUST00000192650.5': [(4491718, 4492668, 2), (4493771, 4493863, 1), (4495135, 4495155, 0)],
                'ENSMUST00000027036.10': [(4807913, 4807982, 0), (4808454, 4808486, 0), (4828583, 4828649, 1),
                                          (4830267, 4830315, 1), (4832310, 4832381, 1), (4837000, 4837074, 2)],
                'ENSMUST00000150971.7': [(4807913, 4807982, 0), (4808454, 4808486, 0), (4828583, 4828649, 1),
                                         (4830267, 4830315, 1), (4832310, 4832381, 1)],
                'ENSMUST00000119612.8': [(4807913, 4807982, 0), (4808454, 4808486, 0), (4828583, 4828649, 1),
                                         (4830267, 4830315, 1), (4832310, 4832368, 1), (4841106, 4841109, 0)],
                }

    def test_gene_absent(self, gene, gene_to_transcript, transcript_to_cds):
        gene.name = 'doesnt exist'
        assert not filter.junction_is_annotated(gene, gene_to_transcript, transcript_to_cds)

    def test_one_junction_present(self, gene, gene_to_transcript, transcript_to_cds):
        # the first transcript of our gene (ENSMUST00000027035.95) has one junction: ['4492668:4493099'];
        # the second transcript of our gene (ENSMUST00000195555.1) has no junction:
        # the third transcript of our gene (ENSMUST00000192650.5) has two junctions: ['4492668:4493771', '4493863:4495135'];
        expected = set(['4492668:4493771', '4493863:4495135', '4492668:4493099'])
        assert expected == filter.junction_is_annotated(gene, gene_to_transcript, transcript_to_cds)


def test_junction_tuple_is_annotated():
    junction_flag = np.zeros((4, 4))
    junction_flag[1, 2] = 1
    junction_flag[2, 3] = 1
    vertex_id_tuple = (1, 2, 3)
    assert (True, True) == filter.junction_tuple_is_annotated(junction_flag, vertex_id_tuple)
    vertex_id_tuple = (0, 2, 3)
    assert (False, True) == filter.junction_tuple_is_annotated(junction_flag, vertex_id_tuple)
    vertex_id_tuple = (0, 1, 3)
    assert (False, False) == filter.junction_tuple_is_annotated(junction_flag, vertex_id_tuple)
    vertex_id_tuple = (1, 2)
    assert (True,) == filter.junction_tuple_is_annotated(junction_flag, vertex_id_tuple)


class TestIsIntronInJunctionList:
    """ Tests for :meth:`filter.is_intron_in_junction_list` """

    def test_empty(self):
        splice_graph = Splicegraph()
        splice_graph.vertices = np.array([[100, 300, 500, 700], [105, 305, 505, 705]], dtype='int')
        vertex_ids = [0, 1, 2, 3]
        assert np.isnan(filter.is_intron_in_junction_list(splice_graph, vertex_ids, '+', None))

    def test_has_nans(self):
        splice_graph = Splicegraph()
        vertex_ids = [np.nan, 1, 2, 3]
        assert np.isnan(filter.is_intron_in_junction_list(splice_graph, vertex_ids, '+', None))

    def test_presence_absence(self):
        splice_graph = Splicegraph()
        splice_graph.vertices = np.array([[100, 300, 500, 700], [105, 305, 505, 705]], dtype='int')
        vertex_ids = [1, 2]
        assert 0 == filter.is_intron_in_junction_list(splice_graph, vertex_ids, '+', [])
        assert 0 == filter.is_intron_in_junction_list(splice_graph, vertex_ids, '+', ['100:200:-', '300:400:+'])
        assert 1 == filter.is_intron_in_junction_list(splice_graph, vertex_ids, '+', ['305:500:+'])
        assert 1 == filter.is_intron_in_junction_list(splice_graph, vertex_ids, '-', ['505:300:-'])
        assert 0 == filter.is_intron_in_junction_list(splice_graph, vertex_ids, '-', ['305:500:+'])


class TestAddKmers:
    """ Tests for :meth:`filter.add_kmers` """

    def test_empty(self):
        kmer_set = set()
        filter.add_kmers(kmer_set, [])
        assert len(kmer_set) == 0

    def test_three_kmers(self):
        kmer_set = set()
        kmers = [OutputKmer(kmer='MSSPDAGYA', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan, junction_annotated=True, reading_frame_annotated=False),
                 OutputKmer(kmer='SSPDAGYAS', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan, junction_annotated=True, reading_frame_annotated=False),
                 OutputKmer(kmer='SPDAGYASD', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan, junction_annotated=True, reading_frame_annotated=False)]

        filter.add_kmers(kmer_set, kmers)
        assert len(kmer_set) == 3
        assert kmer_set == set(['MSSPDAGYA', 'SSPDAGYAS', 'SPDAGYASD'])

    def test_duplicates(self):
        kmer_set = set(['MSSPDAGYA', 'SSPDAGYAS'])
        kmers = [OutputKmer(kmer='MSSPDAGYA', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan, junction_annotated=True, reading_frame_annotated=False),
                 OutputKmer(kmer='SSPDAGYAS', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan, junction_annotated=True, reading_frame_annotated=False),
                 OutputKmer(kmer='SPDAGYASD', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan, junction_annotated=True, reading_frame_annotated=False)]

        filter.add_kmers(kmer_set, kmers)
        assert len(kmer_set) == 3
        assert kmer_set == set(['MSSPDAGYA', 'SSPDAGYAS', 'SPDAGYASD'])


