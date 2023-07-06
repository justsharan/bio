from typing import Generator

NUCLEOTIDE_PAIRS = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

# Parse RNA codon table
CODON_TABLE = {}
with open('RNA_codon_table.txt', 'r') as file:
    for line in file:
        items = line.strip().split()
        CODON_TABLE[items[0]] = items[1] if len(items) == 2 else ''

def window(sequence: str, l: int) -> Generator[str, None, None]:
    """Create a sliding window for `sequence` of length `l`"""
    for i in range(len(sequence)-l+1):
        yield sequence[i:i+l]

def chunks(sequence: str, n: int = 3) -> Generator[str, None, None]:
    """Yields chunks of a sequence at a time"""
    for i in range(0, len(sequence), n):
        yield sequence[i:i+n]

def complement(sequence: str) -> str:
    """Return the complement of a DNA sequence"""
    return ''.join(NUCLEOTIDE_PAIRS[nucleotide] for nucleotide in sequence)

def reverse_comp(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence"""
    return ''.join(NUCLEOTIDE_PAIRS[nucleotide] for nucleotide in reversed(sequence))

def transcribe(sequence: str) -> str:
    """Return the RNA sequence for a given DNA sequence"""
    return sequence.replace('T', 'U')

def translate(sequence: str) -> str:
    """Translate an RNA sequence to a peptide sequence using the standard codon table"""
    return ''.join(CODON_TABLE[codon] for codon in chunks(sequence))

def DNA_to_peptide(sequence: str) -> str:
    """Transcribe and translate a DNA sequence to a peptide"""
    return translate(transcribe(sequence))
