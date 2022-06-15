import numpy as np

from spladder.classes.gene import Gene
from spladder.classes.splicegraph import Splicegraph

from immunopepper.io_ import list_to_tuple, namedtuple_to_str
from immunopepper.namedtuples import Coord, OutputKmer, VertexPair


def _collect_remove_ids(exon_junction_dict):
    """ For all exon pairs around the same intron, keep the longest one.

    :param exon_junction_dict: dictionary mapping (read_frame, junction_start, junction_end) to to the list of
        overlapping exons, i.e. list[tuple(output_idx, pos_start, pos_end)], as created by
        :meth:`_get_exon_junction_dict`
    :return: a list of junction ids to remove
    :rtype: list[str]
    """

    remove_id_list = []
    removed = set()
    for exon_pair_list in exon_junction_dict.values():
        exon_count = len(exon_pair_list)
        if exon_count < 2:
            continue
        for i, exon_pair in enumerate(exon_pair_list):
            if exon_pair[0] in removed:
                continue
            for j in range(i + 1, exon_count):
                if exon_pair_list[j][0] in removed:
                    continue
                i_pos1 = exon_pair[1]
                i_pos2 = exon_pair[2]
                j_pos1 = exon_pair_list[j][1]
                j_pos2 = exon_pair_list[j][2]
                if i_pos1 >= j_pos1 and i_pos2 <= j_pos2:
                    remove_id_list.append(exon_pair[0])
                    removed.add(exon_pair[0])
                    break
                if j_pos1 >= i_pos1 and j_pos2 <= i_pos2:
                    remove_id_list.append(exon_pair_list[j][0])
                    removed.add(exon_pair_list[j][0])

    return remove_id_list


def _get_exon_junction_dict(vertex_pairs: list[VertexPair], strand: str):
    """ Creates an auxiliary dictionary used for filtering. Assigns to each junction all overlapping exon pairs.

    :param vertex_pairs: list of VertexPairs instances used to create the dictionary. From each VertexPair a key and a
        value is generated. The key is reading frame and intron coordinates, the value is the vertex index and exon
        coordinates
    :param strand: '+' or '-' denoting direct or reverse strand
    :return: a dictionary mapping tuples of (read_frame, junction_start, junction_stop) to the list of overlapping
        exons, i.e. list[tuple(output_idx, pos_start, pos_end)]
    :rtype dict:
    """
    exon_dict = {}
    for metadata in vertex_pairs:
        # TODO: need to come up with a new way to index the exon dict
        coord = metadata.modified_exons_coord
        idx = metadata.output_id
        read_frame = metadata.read_frame.read_phase
        # compute intron (i1, i2) and other (o1, o2) coordinates, then add them to the dictionary
        if strand == '+':
            i1, i2, o1, o2 = coord.stop_v1, coord.start_v2, coord.start_v1, coord.stop_v2
        else:
            i1, i2, o1, o2 = coord.stop_v2, coord.start_v1, coord.start_v2, coord.stop_v1
        key = (read_frame, i1, i2)
        if key in exon_dict:
            exon_dict[key].append((idx, o1, o2))
        else:
            exon_dict[key] = [(idx, o1, o2)]

    return exon_dict


def filter_redundant_junctions(vertex_pairs: list[VertexPair], strand: str):
    """ Given a list of exon junctions, remove the ones that redundantly cover a junction, keeping only the longest one

    :param vertex_pairs: list of :class:`VertexPair`s with exon coordinates
    :param strand: strand of the gene, can be '+' or '-'
    :return filtered_vertex_pairs: list of non-redundant :class:`VertexPair`s from vertex_pairs
    :rtype list[VertexPair]:
    """

    exon_dict = _get_exon_junction_dict(vertex_pairs, strand)
    remove_id_list = _collect_remove_ids(exon_dict)
    return list(filter(lambda m: m.output_id not in remove_id_list, vertex_pairs))


def junction_is_annotated(gene: Gene, gene_to_transcript_table: dict[str: list[str]],
                          transcript_to_cds_table: dict[str: list[tuple[int, int, int]]]):
    """ Indicates whether a junction also appears in any transcript given by .gtf file

    :param gene: :class:`Gene` instance created by Spladder
    :param gene_to_transcript_table: maps gene to its transcripts, e.g.
        'ENSMUSG00000025902.13' -> ['ENSMUST00000027035.9', 'ENSMUST00000195555.1']
    :param transcript_to_cds_table: maps a transcript to a list of coding DNA sequencing (cds), e.g.
        'ENSMUST00000027035.9' -> [(4491718, 4492668, 2), (4493099, 4493406, 0)]
    :return junction_flag: 2D array matrix with shape len(edges) x len(edges).
        True indicates the exon pair also appears in gtf file while False indicates that the
         exon pair not exist or it does not show in the given .gtf file.
    """
    vertices = gene.splicegraph.vertices
    edges = gene.splicegraph.edges
    junction_flag = np.zeros(edges.shape, dtype='bool')

    if gene.name not in gene_to_transcript_table:
        return junction_flag

    ts_list = gene_to_transcript_table[gene.name]
    for ts in ts_list:
        if ts not in transcript_to_cds_table:
            continue

        curr_ts = transcript_to_cds_table[ts]
        if len(curr_ts) < 2:
            continue
        curr_ts = np.array(curr_ts)[:, [0, 1]]
        # create a list of pairs that join the end of a transcript to the beginning of the next transcript
        transcript_junctions = [':'.join(x) for x in
                                curr_ts.ravel()[1:-1].reshape(curr_ts.shape[0] - 1, 2).astype('str')]

        for x in np.array(np.where(np.triu(edges))).T:
            if f'{vertices[1, x[0]]}:{vertices[0, x[1]]}' in transcript_junctions:
                junction_flag[x[0], x[1]] = 1
                junction_flag[x[1], x[0]] = 1

    return junction_flag


def junction_tuple_is_annotated(junction_flag: np.ndarray, vertex_ids: list[int]):
    """Check if any junction in the given list vertices is flagged.

    :param junction_flag: 2-d array with shape (number_of_vertices, number_of_vertices). junction_flag[i,j] == 1 iff
       vertices i and j are connected
    :param vertex_ids: list of vertex ids defining junction pairs
    :return:
       np.nan if vertex_ids contains a np.nan,
       otherwise returns a tuple containing 0 or 1 flag for each of the junctions. Example (1,) or (1,0) or (1,1)
    """
    if np.nan in vertex_ids:
        return (np.nan)

    junction_annotated = [1 if junction_flag[vertex_ids[i], vertex_ids[i + 1]] else 0
                          for i in range(len(vertex_ids) - 1)]
    return tuple(junction_annotated)


def is_intron_in_junction_list(splicegraph: Splicegraph, vertex_ids: list[int], strand: str,
                               junction_list: list[str]):
    """
    Check if the intron defined by vertex_ids is in the user provided list of junctions
    :param splicegraph: Spladder splice graph
    :param vertex_ids: list vertex coordinates in splicegraph defining the intron (typically 2 coordinates)
    :param strand: '+' for direct, '-' for reverse strand
    :param junction_list: list of junctions to check the intron against
    :return: np.nan if junction_list is None or if vertex_ids contains a np.nan
             1 if the intron is in junction_list, 0 if not
    """

    if np.nan in vertex_ids or junction_list is None:
        return np.nan

    vertex_coord_list = splicegraph.vertices
    if strand == '-':
        vertex_ids = vertex_ids[::-1]  # in negative strand case, the vertex id list is [5,3] instead of [3,5]
    for i in range(len(vertex_ids) - 1):
        junction_comb = f'{vertex_coord_list[1, vertex_ids[i]]}:{vertex_coord_list[0, vertex_ids[i + 1]]}:{strand}'
        if junction_comb in junction_list:
            return 1
    return 0


# TODO(dd): this function smells. It takes nicely structured data and massages it into an undecipherable jumble. Why?
def add_kmer_properties(foreground_dict, output_kmers: list[OutputKmer]):
    """
    Update foreground_dict by collapsing the metadata in output_kmers.

    :param foreground_dict: a dictionary where keys are peptides, values are lists of property sets containing the
        collapsed metadata for the given peptide; for example::
        {'MSPHGYKLA' : [{'ENSMUSG00000025902.13'}, {51.43}, {True}, {1.0}]}
    :param output_kmers: list of :class:`OutputKmer` used to update the dictionary
    """
    for output_kmer in output_kmers:
        # Prepare metadata
        ord_dict_metadata = output_kmer._asdict()
        del ord_dict_metadata['kmer']

        # aggregate metadata of unique kmers
        if output_kmer.kmer not in foreground_dict:
            add_novel_kmer = []
            for i, j in enumerate(ord_dict_metadata.values()):
                if i == 0:
                    add_novel_kmer.append({j.split(':')[0]})
                else:
                    add_novel_kmer.append({j})
            foreground_dict[output_kmer.kmer] = add_novel_kmer
        else:
            for idx, value in enumerate(ord_dict_metadata.values()):
                if idx == 0:
                    foreground_dict[output_kmer.kmer][idx].add(value.split(':')[0])
                else:
                    foreground_dict[output_kmer.kmer][idx].add(value)


def add_kmers(kmer_set: set[str], output_kmers: list[OutputKmer]):
    """ Add the kmer from each output_kmers to kmer_set. """
    for output_kmer in output_kmers:
        kmer_set.add(output_kmer.kmer)


# TODO(dd): this function also smells. It takes nicely structured data and massages it into an undecipherable jumble
def add_peptide_properties(peptide_properties: dict[str, list[set]], namedtuple_list: list,
                           skip_expr=False):
    """
    Updates peptide_properties using namedtuple_list. The namedtuple_list is converted to a dictionary, and each element
    from this dictionary is used to update the corresponding values for the given peptide. Lists and nested named tuples
    within namedtuple_list are converted to tuples prior to addition to the set.

    :param peptide_properties: maps peptides to a set of properties for that peptide (e.g. 'output_id`, 'chr',
          'gene_name', etc). The properties are extracted from the namedtuple_list.
    :param namedtuple_list: list of `class:OutputBackground` or of `class:OutputMetadata`
    :param skip_expr: if True, the 'junction_expr' and 'segment_expr' keys are removed from the result
    """
    for namedtuple_peptide in namedtuple_list:
        ord_dict_metadata = {}
        for k, v in namedtuple_peptide._asdict().items():
            ord_dict_metadata[k] = namedtuple_to_str(v, None, ';') if isinstance(v, Coord) else list_to_tuple(v);

        del ord_dict_metadata['peptide']
        if skip_expr:
            del ord_dict_metadata['junction_expr']
            del ord_dict_metadata['segment_expr']

        # aggregate metadata of unique peptides
        if namedtuple_peptide.peptide not in peptide_properties:
            peptide_properties[namedtuple_peptide.peptide] = [{i} for i in ord_dict_metadata.values()]
        else:
            for idx, value in enumerate(ord_dict_metadata.values()):
                peptide_properties[namedtuple_peptide.peptide][idx].add(value)
