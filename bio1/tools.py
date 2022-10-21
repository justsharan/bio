import numpy as np

def print_motifs(motifs: np.ndarray, sep='\n'):
  final_str = ''
  for motif in motifs:
    final_str += ''.join(motif)
    final_str += sep
  print(final_str)

def window(text: str, l: int):
  for i in range(0, len(text)-l+1):
    yield text[i:i+l]