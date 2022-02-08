"""Contains the functionality used for filtering"""

import numpy as np

from .io_ import convert_to_str_Coord_namedtuple
from .io_ import list_to_tuple


def _collect_remove_ids(exon_junction_dict):
    """ For all exon pairs around the same intron, keep the longest one.

    :dict exon_junction_dict: Dictionary mapping (read_frame, pos_mid_1, pos_mid_2, variant_combination)
                               -> List[Tuple(output_idx, pos_start, pos_end,)]
    :return: a list of junction ids to remove
    :rtype: list
    """

    # TODO: [Warning!] if two output lines have identical coordinates and readframe both of them will be removed
    remove_id_list = []
    removed = set()
    for exon_pair_list in exon_junction_dict.values():
        size = len(exon_pair_list)
        if size < 2:
            continue
        for i, exon_pair in enumerate(exon_pair_list):
            if exon_pair[0] in removed:
                continue
            for j in range(i + 1, size):
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


def _get_exon_junction_dict(metadata_list, strand):
    """ Creates an auxiliary dictionary used for filtering. Assigns to each junction all overlapping exon pairs.

    :param metadata_list: List(Output_metadata).
    :return exon_dict: dict. (read_frame, pos_mid_1, pos_mid_2, variant_combination)
        -> List[Tuple(output_idx, pos_start, pos_end,)]
    :rtype dict:
    """
    exon_dict = {}
    for metadata in metadata_list:
        ## TODO: need to come up with a new way to index the exon dict
        coord = metadata.modified_exons_coord
        idx = metadata.output_id
        read_frame = metadata.read_frame.read_phase
        ### figure out intron vs other coordinates
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


def get_filtered_metadata_list(metadata_list, strand):
    """ Given a list of exon junctions, remove the ones that redundantly cover a junction

    :param metadata_list: List(Output_metadata),
    :param strand: strand of the gene

    :return filtered_meetadata_list: List of metadata objects remaining after filtering
    :rtype list:
    """

    exon_dict = _get_exon_junction_dict(metadata_list, strand)
    remove_id_list = _collect_remove_ids(exon_dict)
    return list(filter(lambda m: m.output_id not in remove_id_list, metadata_list))


def junction_is_annotated(gene, gene_to_transcript_table, transcript_to_cds_table):
    """ Indicates whether a junction also appears in any transcript given by .gtf file

    :param gene: Object. Created by SplAdder.
    :param gene_to_transcript_table: Dict. "Gene" -> "transcript"
    :param transcript_to_cds_table: Dict. "transcript" -> "cds"
    :return junction_flag: 2D array matrix with shape len(edges)xlen(edges).
        True indicates the exon pair also appear in gtf file while False indicates the
         exon pair not exist or it does not show in the given .gtf file.
    """
    vertices = gene.splicegraph.vertices
    edges = gene.splicegraph.edges
    junction_flag = np.zeros(edges.shape, dtype='bool')

    if gene.name not in gene_to_transcript_table:
        return junction_flag

    ts_list = gene_to_transcript_table[gene.name]
    for ts in ts_list:
        if not ts in transcript_to_cds_table:
            continue

        curr_ts = transcript_to_cds_table[ts]
        if len(curr_ts) < 2:
            continue
        curr_ts = np.array(curr_ts)[:, [0, 1]]
        transcript_junctions = [':'.join(x) for x in
                                curr_ts.ravel()[1:-1].reshape(curr_ts.shape[0] - 1, 2).astype('str')]

        for x in np.array(np.where(np.triu(edges))).T:
            if '%i:%i' % (vertices[1, x[0]], vertices[0, x[1]]) in transcript_junctions:
                junction_flag[x[0], x[1]] = 1
                junction_flag[x[1], x[0]] = 1

    return junction_flag


def junction_tuple_is_annotated(junction_flag, vertex_id_tuple):
    """Returns 1 if any junction in the given tuple of vertexes is flagged and 0 otherwise.

    :param junction_flag: 2-d array with shape (number_of_vertices, number_of_vertices). junction_flag(i,j) is 1 iff
       vertices i and j are connected
    :param vertex_id_tuple: tuple of vertex ids defining junction pairs
    :return:
       np.nan if vertix_id_tuple contains a np.nan
       1 if any junction in param:`vertex_id_tuple` is connected in param:`junction_flag`, 0 otherwise
    """
    if np.nan in vertex_id_tuple:
        return np.nan
    for i in range(len(vertex_id_tuple) - 1):
        if junction_flag[vertex_id_tuple[i], vertex_id_tuple[i + 1]]:
            return 1
    return 0


def junction_is_in_given_list(splicegraph, vertex_id_list, strand, junction_list):
    """Check if the intron is in the user provided list of junctions"""

    if np.nan in vertex_id_list or junction_list is None:
        return np.nan

    vertex_coord_list = splicegraph.vertices
    if strand == '-':
        vertex_id_list = vertex_id_list[::-1]  # in negative strand case, the vertex id list is [5,3] instead of [3,5]
    for i in range(len(vertex_id_list) - 1):
        junc_comb = ':'.join(
            [str(vertex_coord_list[1, vertex_id_list[i]]), str(vertex_coord_list[0, vertex_id_list[i + 1]]), strand])
        if junc_comb in junction_list:
            return 1
    return 0


def add_dict_kmer_forgrd(foregrd_dict, _namedtuple_list):
    """...
        From a namedtuple
        Updates the following dictionnary:
        keys = unique peptides,
        values = collapsed metadata, each metadata item is converted to a set. {'output_id' : set(), 'chr' : set()}

        Parameters
        ----------
        foregrd_dict:
        _namedtuple_list:

        """
    for _namedtuple_kmer in _namedtuple_list:
        ### Prepare metadata
        ord_dict_metadata = _namedtuple_kmer._asdict()
        del ord_dict_metadata['kmer']

        ### aggregate metadata of unique kmers
        add_novel_kmer = []
        if _namedtuple_kmer.kmer not in foregrd_dict:
            for i, j in enumerate(ord_dict_metadata.values()):
                if i == 0:
                    add_novel_kmer.append({j.split(':')[0]})
                else:
                    add_novel_kmer.append({j})
            foregrd_dict[_namedtuple_kmer.kmer] = add_novel_kmer
        else:
            for id, value in enumerate(ord_dict_metadata.values()):
                if id == 0:
                    foregrd_dict[_namedtuple_kmer.kmer][id].add(value.split(':')[0])
                else:
                    foregrd_dict[_namedtuple_kmer.kmer][id].add(value)


def add_set_kmer_back(background_dict, namedtuple_list):
    """
        Adds the kmer field from namedtuple_list to background_dict, if not already present.
        :param background_dict:
        :param namedtuple_list:
    """
    for _namedtuple_kmer in namedtuple_list:
        if _namedtuple_kmer.kmer in background_dict:
            continue
        background_dict.add(_namedtuple_kmer.kmer)


def add_dict_peptide(dict_peptides, namedtuple_list, skip_expr=False):
    """
    Updates dict_peptides using namedtuple_list.
        keys = unique peptides,
        values = collapsed metadata, each metadata item is converted to a set. {'output_id' : set(), 'chr' : set()}
        Lists and nested namedtuples are converted to tuples prior to addition to the set.

        :param dict_peptides:
        :param namedtuple_list:
        :param skip_expr: if True, the 'junction_expr' and 'segment_expr' keys are removed from the result

    """
    for namedtuple_peptide in namedtuple_list:
        ord_dict_metadata = dict(namedtuple_peptide._asdict())
        ord_dict_metadata = {k: list_to_tuple(v) for k, v in ord_dict_metadata.items()}
        ord_dict_metadata = {k: (convert_to_str_Coord_namedtuple(v, ';')) for k, v in ord_dict_metadata.items()}
        del ord_dict_metadata['peptide']
        if skip_expr:
            del ord_dict_metadata['junction_expr']
            del ord_dict_metadata['segment_expr']

        ### aggregate metadata of unique peptides
        if namedtuple_peptide.peptide not in dict_peptides:
            dict_peptides[namedtuple_peptide.peptide] = [{i} for i in ord_dict_metadata.values()]
        else:
            for id, value in enumerate(ord_dict_metadata.values()):
                dict_peptides[namedtuple_peptide.peptide][id].add(value)
