
seq = "A T T G C C C C G A A T"

alphabet = set()

for c in seq:
    if c != ' ':
        alphabet.add(c)

print("Alphabet of the sequence:", " ".join(sorted(alphabet)))
