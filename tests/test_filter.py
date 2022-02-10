"""
Unit tests for imunopepper/filter.py
"""
import numpy as np
import pytest
from spladder.classes.gene import Gene
from spladder.classes.splicegraph import Splicegraph

import immunopepper.filter as filter
from immunopepper.namedtuples import Coord, OutputKmer, OutputMetadata, ReadingFrameTuple, VertexPair


class TestFilterRedundantJunctions:
    """ Tests for :meth:`filter.filter_redundant_junctions` """
    def test_no_redundant(self):
        read_frame = ReadingFrameTuple(cds_left_modi=4495135, cds_right_modi=4495155, read_phase=2)
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
        read_frame = ReadingFrameTuple(cds_left_modi=4495135, cds_right_modi=4495155, read_phase=2)
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
        read_frame = ReadingFrameTuple(cds_left_modi=4495135, cds_right_modi=4495155, read_phase=2)
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
        assert not np.any(filter.junction_is_annotated(gene, gene_to_transcript, transcript_to_cds))

    def test_gene_present_empty_graph(self, gene, gene_to_transcript, transcript_to_cds):
        assert not np.any(filter.junction_is_annotated(gene, gene_to_transcript, transcript_to_cds))

    def test_one_junction_present(self, gene, gene_to_transcript, transcript_to_cds):
        # the third transcript of our gene (ENSMUST00000192650.5) has two junctions:
        # ['4492668:4493771', '4493863:4495135']; let's create a graph that has one of these junctions
        gene.splicegraph.vertices = np.array([[4491718, 4493771, 4495135], [123456, 4493863, 4495155]], dtype='int')
        gene.splicegraph.edges = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype='int')
        expected = np.array([[False, False, False], [False, False, True], [False, True, False]])
        assert np.array_equal(expected, filter.junction_is_annotated(gene, gene_to_transcript, transcript_to_cds))

    def test_two_junctions_present(self, gene, gene_to_transcript, transcript_to_cds):
        # the third transcript of our gene (ENSMUST00000192650.5) has two junctions:
        # ['4492668:4493771', '4493863:4495135']; let's create a graph that has *both* of these junctions
        gene.splicegraph.vertices = np.array([[4491718, 4493771, 4495135], [4492668, 4493863, 4495155]], dtype='int')
        gene.splicegraph.edges = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype='int')
        expected = np.array([[False, True, False], [True, False, True], [False, True, False]])
        assert np.array_equal(expected, filter.junction_is_annotated(gene, gene_to_transcript, transcript_to_cds))


def test_junction_tuple_is_annotated():
    junction_flag = np.zeros((4, 4))
    junction_flag[1, 2] = 1
    junction_flag[2, 3] = 1
    vertex_id_tuple = (1, 2, 3)
    assert 1 == filter.junction_tuple_is_annotated(junction_flag, vertex_id_tuple)
    vertex_id_tuple = (0, 2, 3)
    assert 1 == filter.junction_tuple_is_annotated(junction_flag, vertex_id_tuple)
    vertex_id_tuple = (0, 1, 3)
    assert 0 == filter.junction_tuple_is_annotated(junction_flag, vertex_id_tuple)
    vertex_id_tuple = (1, 2)
    assert 1 == filter.junction_tuple_is_annotated(junction_flag, vertex_id_tuple)


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


class TestAddKmerProperties:
    """ Tests for :meth:`filter.add_kmer_properties` """

    def test_add(self):
        output_kmers = [
            OutputKmer(kmer='MSPHGYKLA', id='ENSMUSG00000025902.13:32_14:0:4495155:2-exons', segment_expr=51.43,
                       is_cross_junction=True, junction_expr=1.0),
            OutputKmer(kmer='SPHGYKLAS', id='ENSMUSG00000025902.13:32_14:0:4495155:2-exons', segment_expr=57.45,
                       is_cross_junction=True, junction_expr=1.0),
            OutputKmer(kmer='PHGYKLASD', id='ENSMUSG00000025902.13:32_14:0:4495155:2-exons', segment_expr=63.48,
                       is_cross_junction=True, junction_expr=1.0),
            OutputKmer(kmer='HGYKLASDP', id='ENSMUSG00000025902.13:32_14:0:4495155:2-exons', segment_expr=69.5,
                       is_cross_junction=True, junction_expr=1.0)]
        kmer_dict = {}
        filter.add_kmer_properties(kmer_dict, output_kmers)
        assert len(kmer_dict) == 4
        assert kmer_dict == {'MSPHGYKLA': [{'ENSMUSG00000025902.13'}, {51.43}, {True}, {1.0}],
                             'SPHGYKLAS': [{'ENSMUSG00000025902.13'}, {57.45}, {True}, {1.0}],
                             'PHGYKLASD': [{'ENSMUSG00000025902.13'}, {63.48}, {True}, {1.0}],
                             'HGYKLASDP': [{'ENSMUSG00000025902.13'}, {69.5}, {True}, {1.0}]}

    def test_duplicate_kmer(self):
        output_kmers = [
            OutputKmer(kmer='MSPHGYKLA', id='ENSMUSG00000025902.13:32_14:0:4495155:2-exons', segment_expr=51.43,
                       is_cross_junction=False, junction_expr=1.0),
            OutputKmer(kmer='MSPHGYKLA', id='ENSMUSG00000025902.13:32_14:0:4495155:2-exons', segment_expr=12.34,
                       is_cross_junction=True, junction_expr=1.0),
            OutputKmer(kmer='SPHGYKLAS', id='ENSMUSG00000025902.13:32_14:0:4495155:2-exons', segment_expr=57.45,
                       is_cross_junction=True, junction_expr=1.0),
            OutputKmer(kmer='PHGYKLASD', id='ENSMUSG00000025902.13:32_14:0:4495155:2-exons', segment_expr=63.48,
                       is_cross_junction=True, junction_expr=1.0),
            OutputKmer(kmer='HGYKLASDP', id='ENSMUSG00000025902.13:32_14:0:4495155:2-exons', segment_expr=69.5,
                       is_cross_junction=True, junction_expr=1.0)]
        kmer_dict = {}
        filter.add_kmer_properties(kmer_dict, output_kmers)
        assert len(kmer_dict) == 4
        assert kmer_dict == {'MSPHGYKLA': [{'ENSMUSG00000025902.13'}, {51.43, 12.34}, {True, False}, {1.0}],
                             'SPHGYKLAS': [{'ENSMUSG00000025902.13'}, {57.45}, {True}, {1.0}],
                             'PHGYKLASD': [{'ENSMUSG00000025902.13'}, {63.48}, {True}, {1.0}],
                             'HGYKLASDP': [{'ENSMUSG00000025902.13'}, {69.5}, {True}, {1.0}]}


class TestAddKmers:
    """ Tests for :meth:`filter.add_kmers` """

    def test_empty(self):
        kmer_set = set()
        filter.add_kmers(kmer_set, [])
        assert len(kmer_set) == 0

    def test_three_kmers(self):
        kmer_set = set()
        kmers = [OutputKmer(kmer='MSSPDAGYA', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan),
                 OutputKmer(kmer='SSPDAGYAS', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan),
                 OutputKmer(kmer='SPDAGYASD', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan)]

        filter.add_kmers(kmer_set, kmers)
        assert len(kmer_set) == 3
        assert kmer_set == set(['MSSPDAGYA', 'SSPDAGYAS', 'SPDAGYASD'])

    def test_duplicates(self):
        kmer_set = set(['MSSPDAGYA', 'SSPDAGYAS'])
        kmers = [OutputKmer(kmer='MSSPDAGYA', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan),
                 OutputKmer(kmer='SSPDAGYAS', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan),
                 OutputKmer(kmer='SPDAGYASD', id='ENSMUST00000027035.9', segment_expr=97.29, is_cross_junction=False,
                            junction_expr=np.nan)]

        filter.add_kmers(kmer_set, kmers)
        assert len(kmer_set) == 3
        assert kmer_set == set(['MSSPDAGYA', 'SSPDAGYAS', 'SPDAGYASD'])


class TestAddPeptideProperties:
    """ Tests for :meth:`filter.add_peptide_properties` """

    def test_empty(self):
        peptide_properties = {}
        output_metadata = []
        filter.add_peptide_properties(peptide_properties, output_metadata)
        assert len(peptide_properties) == 0

    def test_add_one(self):
        peptide_properties = {}
        output_metadata = [
            OutputMetadata(peptide='MSPHGYK', output_id='ENSMUSG00000025902.13:23_4:0:4495155:2-exons', read_frame=1,
                           gene_name='ENSMUSG00000025902.13', gene_chr='chr1', gene_strand='-', mutation_mode='ref',
                           junction_annotated=1, has_stop_codon=1, is_in_junction_list=np.nan, is_isolated=1,
                           variant_comb=np.nan,
                           variant_seg_expr=np.nan,
                           modified_exons_coord=Coord(start_v1=4493771, stop_v1=4495155, start_v2=4491712,
                                                      stop_v2=4492668,
                                                      start_v3=None, stop_v3=None),
                           original_exons_coord=Coord(start_v1=4493771, stop_v1=4496413, start_v2=4491712,
                                                      stop_v2=4492668,
                                                      start_v3=None, stop_v3=None), vertex_idx=[23, 4],
                           junction_expr=np.nan,
                           segment_expr=99, kmer_type='2-exons')]
        filter.add_peptide_properties(peptide_properties, output_metadata)
        assert len(peptide_properties) == 1
        assert list(peptide_properties.keys()) == ['MSPHGYK']
        assert peptide_properties['MSPHGYK'] == [{'ENSMUSG00000025902.13:23_4:0:4495155:2-exons'}, {1},
                                                 {'ENSMUSG00000025902.13'}, {'chr1'}, {'-'}, {'ref'}, {1}, {1},
                                                 {np.nan},
                                                 {1},
                                                 {np.nan}, {np.nan}, {'4493771;4495155;4491712;4492668;None;None'},
                                                 {'4493771;4496413;4491712;4492668;None;None'}, {(23, 4)}, {np.nan},
                                                 {99},
                                                 {'2-exons'}]

    def test_add_two(self):
        peptide_properties = {}
        output_metadata1 = [
            OutputMetadata(peptide='MSPHGYK', output_id='ENSMUSG00000025902.13:23_4:0:4495155:2-exons', read_frame=1,
                           gene_name='ENSMUSG00000025902.13', gene_chr='chr1', gene_strand='-', mutation_mode='ref',
                           junction_annotated=1, has_stop_codon=1, is_in_junction_list=np.nan, is_isolated=1,
                           variant_comb=np.nan,
                           variant_seg_expr=np.nan,
                           modified_exons_coord=Coord(start_v1=4493771, stop_v1=4495155, start_v2=4491712,
                                                      stop_v2=4492668,
                                                      start_v3=None, stop_v3=None),
                           original_exons_coord=Coord(start_v1=4493771, stop_v1=4496413, start_v2=4491712,
                                                      stop_v2=4492668,
                                                      start_v3=None, stop_v3=None), vertex_idx=[23, 4],
                           junction_expr=np.nan,
                           segment_expr=99, kmer_type='2-exons')]
        output_metadata2 = [
            OutputMetadata(peptide='MSPHGYK', output_id='ENSMUSG00000025902.13:23_4:0:4495155:2-exons2', read_frame=1,
                           gene_name='ENSMUSG00000025902.13', gene_chr='chr2', gene_strand='+', mutation_mode='ref',
                           junction_annotated=1, has_stop_codon=1, is_in_junction_list=np.nan, is_isolated=1,
                           variant_comb=np.nan,
                           variant_seg_expr=np.nan,
                           modified_exons_coord=Coord(start_v1=4493771, stop_v1=4495155, start_v2=4491712,
                                                      stop_v2=4492668,
                                                      start_v3=None, stop_v3=None),
                           original_exons_coord=Coord(start_v1=4493771, stop_v1=4496413, start_v2=4491712,
                                                      stop_v2=4492668,
                                                      start_v3=None, stop_v3=None), vertex_idx=[23, 4],
                           junction_expr=np.nan,
                           segment_expr=99, kmer_type='2-exons')]
        filter.add_peptide_properties(peptide_properties, output_metadata1)
        filter.add_peptide_properties(peptide_properties, output_metadata2)
        assert len(peptide_properties) == 1
        assert list(peptide_properties.keys()) == ['MSPHGYK']
        assert peptide_properties['MSPHGYK'] == [
            {'ENSMUSG00000025902.13:23_4:0:4495155:2-exons', 'ENSMUSG00000025902.13:23_4:0:4495155:2-exons2'}, {1},
            {'ENSMUSG00000025902.13'}, {'chr1', 'chr2'}, {'-', '+'}, {'ref'}, {1}, {1}, {np.nan},
            {1},
            {np.nan}, {np.nan}, {'4493771;4495155;4491712;4492668;None;None'},
            {'4493771;4496413;4491712;4492668;None;None'}, {(23, 4)}, {np.nan}, {99},
            {'2-exons'}]
