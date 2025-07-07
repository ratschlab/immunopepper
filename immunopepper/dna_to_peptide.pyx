def dna_to_peptide(dna_str: str, all_read_frames: bool):
    """ Translate a DNA sequence encoding a peptide to the corresponding amino-acid sequence.

    All codons that include an 'N' will be translated to 'X'.

    :param dna_str: dna string to be translated.
    :param all_read_frames: if false, only the first peptide until the stop codon is returned, otherwise a list of
      all translated peptides (for each stop codon) is provided
    :return: a tuple containing the translated peptide and a boolean indicating if a stop codon was found
    """
    cdef dict codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'
    }
    dna_str = dna_str.upper()
    cdef list peptides = []
    cdef bint has_stop_codon = False
    cdef list aa_str = []
    for idx in range(0, (len(dna_str) // 3) * 3, 3):
        codon = dna_str[idx:idx + 3]
        if 'N' in codon:
            aa_str.append('X')
        else:
            amino_acid = codontable[codon]
            if amino_acid == '_':
                has_stop_codon = True
                if not all_read_frames:
                    return [''.join(aa_str)], has_stop_codon
                peptides.append(''.join(aa_str))
                aa_str.clear()
            else:
                aa_str.append(amino_acid)

    peptides.append(''.join(aa_str))
    return peptides, has_stop_codon
