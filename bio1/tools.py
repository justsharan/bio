def window(text: str, l: int):
  for i in range(0, len(text)-l+1):
    yield text[i:i+l]