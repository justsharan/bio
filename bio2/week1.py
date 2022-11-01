def window(text: str, l: int):
  for i in range(0, len(text)-l+1):
    yield text[i:i+l]

def composition(text: str, k: int) -> list[str]:
  return list(window(text, k))

def path_to_genome(genome_path: list[str]) -> str:
  return genome_path[0] + ''.join([pattern[-1] for pattern in genome_path[1:]])

def overlap_graph(patterns: list[str]) -> list[(str, list[str])]:
  final = []
  for pattern in patterns:
    matches = list(filter(lambda seq: seq[:-1] == pattern[1:], patterns))
    if len(matches) > 0:
      final.append((pattern, matches))
  return final

def de_bruijn(k: int, text: str) -> dict:
  de_bruijn_map = {}
  for node in set(window(text, k-1)):
    vals = sorted([kmer[1:] for kmer in window(text, k) if kmer[:-1] == node])
    if len(vals) > 0:
      de_bruijn_map[node] = vals
  return de_bruijn_map

def de_bruijn_from_kmers(patterns: list[str]) -> dict:
  de_bruijn_map = {}
  for item in set(pattern[:-1] for pattern in patterns):
    vals = [kmer[1:] for kmer in patterns if kmer[:-1] == item]
    if len(vals) > 0:
      de_bruijn_map[item] = sorted(vals)
  return de_bruijn_map
