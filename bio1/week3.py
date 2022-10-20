from itertools import product
from tools import window
from week2 import hamming_distance
import numpy as np

def combination(k: int):
    return (''.join(p) for p in product('ATCG', repeat=k))

def motif_enumeration(dna: str, k: int, d: int):
  patterns = set()
  for combo in combination(k):
    if all(any(hamming_distance(combo, pat) <= d for pat in window(string, k)) for string in dna):
      patterns.add(combo)
  return patterns

def count(motifs: np.ndarray) -> np.ndarray:
  res = np.zeros((4, motifs.shape[1]), dtype=np.int32)
  for i, el in enumerate('ACGT'):
    res[i,:] = (motifs==el).sum(axis=0)
  return res

def profile(motifs: np.ndarray, pseudocount=0):
  res = count(motifs) + pseudocount
  return res / res.sum(axis=0)

def score(motifs: np.ndarray) -> int:
  res = count(motifs)
  scores = res.sum(axis=0) - res.max(axis=0)
  return scores.sum()

def entropy(col: np.ndarray) -> float:
  col = col[col != 0]
  return (-col * np.log2(col)).sum()

def d(pattern: str, dna: list[str]):
  hd_sum = 0
  for strand in dna:
    distances = [hamming_distance(kmer, pattern) for kmer in window(strand, len(pattern))]
    hd_sum += min(distances)
  return hd_sum

def median_string(dna: list[str], k: int) -> tuple[str, int]:
  (min_kmer, min_d) = ('', np.inf)
  for kmer in window(''.join(dna), k):
    hd = d(kmer, dna)
    if hd < min_d:
      (min_kmer, min_d) = (kmer, hd)
  return (min_kmer, min_d)

def profile_prob(profile_matrix: np.ndarray, kmer: str) -> float:
  prob = 1
  for i, s in enumerate(kmer):
    prob *= profile_matrix['ACGT'.index(s),i]
  return prob

def profitable_kmer(text: str, k: int, profile_matrix: np.ndarray) -> str:
  profits = {kmer: profile_prob(profile_matrix, kmer) for kmer in window(text, k)}
  return max(profits, key=profits.get)

def greedy_motif_search(dna: list[str], k: int, t: int) -> np.ndarray:
  best_motifs = np.array([list(strand[0:k]) for strand in dna])
  for kmer in window(dna[0], k): # only go through kmers in first strand
    motifs = np.array([list(kmer)])
    for strand in dna[1:t]:
      most_profit_kmer = profitable_kmer(strand, k, profile(motifs, pseudocount=1))
      motifs = np.vstack([motifs, list(most_profit_kmer)])
    if score(motifs) < score(best_motifs):
      best_motifs = motifs
  return best_motifs

def distance_between_pattern_and_strings(pattern: str, dna: list[str]) -> int:
  distance = 0
  for strand in dna:
    hamm_dist = np.inf
    for kmer in window(strand, len(pattern)):
      hamm_dist = min(hamm_dist, hamming_distance(pattern, kmer))
    distance += hamm_dist
  return distance
