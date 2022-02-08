import pytest

import immunopepper.dna_to_peptide as dtp


@pytest.fixture
def empty_string():
    return dtp.dna_to_peptide('', False)


def test_empty_string(empty_string):
    assert empty_string[0] == ['']


def test_one_aminoacid():
    peptides, stop = dtp.dna_to_peptide('ATA', False)
    assert len(peptides) == 1
    assert peptides == ['I']


def test_one_aminoacid_one_extra_base():
    peptides, stop = dtp.dna_to_peptide('ATAC', False)
    assert len(peptides) == 1
    assert peptides == ['I']


def test_one_aminoacid_two_extra_bases():
    peptides, stop = dtp.dna_to_peptide('ATACC', False)
    assert len(peptides) == 1
    assert peptides == ['I']


def test_stop_codon_only():
    peptides, stop = dtp.dna_to_peptide('TAA', False)
    assert len(peptides) == 1
    assert peptides == ['']


def test_one_invalid():
    peptides, stop = dtp.dna_to_peptide('ANA', False)
    assert len(peptides) == 1
    assert peptides == ['X']


def test_one_peptide_lowercase():
    peptides, stop = dtp.dna_to_peptide('aaatggcctcatggttacgaattctatggaacggctgctaacgttatctgccggcgtgct', True)
    assert len(peptides) == 1
    assert peptides == ['KWPHGYEFYGTAANVICRRA']


def test_one_peptide_lowercase_extra_bases():
    peptides, stop = dtp.dna_to_peptide('aaatggcctcatggttacgaattctatggaacggctgctaacgttatctgccggcgtgctCC', False)
    assert len(peptides) == 1
    assert peptides == ['KWPHGYEFYGTAANVICRRA']


def test_one_peptide_uppercase():
    peptides, stop = dtp.dna_to_peptide('aaatggcctcatggttacgaattctatggaacggctgctaacgttatctgccggcgtgct', True)
    assert len(peptides) == 1
    assert peptides == ['KWPHGYEFYGTAANVICRRA']


def test_one_peptide_mixedcase():
    peptides, stop = dtp.dna_to_peptide('aaatggcctcatGGTtACGaattctatggaacggctgctaacgttatctgccggcgtgct', False)
    assert len(peptides) == 1
    assert peptides == ['KWPHGYEFYGTAANVICRRA']


def test_two_peptides():
    peptides, stop = dtp.dna_to_peptide('aaatggcctcatggttacgaattctatggaacggctgctaacgttatctgccggcgtgctTAA'
                                        'aaatggcctcatggttacgaattctatggaacggctgctaacgttatctgccggcgtgct', True)
    assert len(peptides) == 2
    assert peptides == ['KWPHGYEFYGTAANVICRRA', 'KWPHGYEFYGTAANVICRRA']

def test_two_peptides_stop_at_1():
    peptides, stop = dtp.dna_to_peptide('aaatggcctcatggttacgaattctatggaacggctgctaacgttatctgccggcgtgctTAA'
                                        'aaatggcctcatggttacgaattctatggaacggctgctaacgttatctgccggcgtgct', False)
    assert len(peptides) == 1
    assert peptides == ['KWPHGYEFYGTAANVICRRA']
