def translate_dna_to_peptide(dna_str:str,all_read_frames:bool):
    """ Translate a DNA sequence encoding a peptide to amino-acid sequence via RNA.

    If 'N' is included in input dna, 'X' will be outputted since 'N' represents
    uncertainty. Also will output a flag indicating if has stop codon.

    Parameters
    ----------
    dna_str: str or List(str). dna string to be translated.
    all_read_frames: boolean for all reading frames translation versus annotated reading frame propagation. Given as string

    Returns
    -------
    aa_str: translated peptide
    has_stop_codon: Indicator for showing if the input dna contains stop codon

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
    cdef list multiple_pep = []
    cdef bint has_stop_codon = False
    cdef list aa_str = []
    for idx in range(0, len(dna_str), 3):
        codon = dna_str[idx:idx + 3]
        if len(codon) < 3:
            break
        if 'N' in codon:
            aa_str.append('X')
        else:
            if codontable[codon] == '_':
                has_stop_codon = True
                if not all_read_frames:
                    return ''.join(aa_str), has_stop_codon
                else:
                    multiple_pep.append(''.join(aa_str))
                    aa_str.clear()
            else:
                aa_str.append(codontable[codon])

    multiple_pep.append(''.join(aa_str))
    return multiple_pep, has_stop_codon
