from typing import Generator

NUCLEOTIDE_PAIRS = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

def window(sequence: str, l: int) -> Generator[str, None, None]:
    """Create a sliding window for `sequence` of length `l`"""
    for i in range(len(sequence)-l+1):
        yield sequence[i:i+l]

def complement(sequence: str) -> str:
    """Return the complement of a DNA sequence"""
    return ''.join(NUCLEOTIDE_PAIRS[nucleotide] for nucleotide in sequence)

def reverse_comp(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence"""
    return ''.join(NUCLEOTIDE_PAIRS[nucleotide] for nucleotide in reversed(sequence))
