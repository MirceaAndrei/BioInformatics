seq = input("Enter a sequence: ")

for symbol in set(seq):
    freq = seq.count(symbol) / len(seq)
    print(symbol, ":", round(freq, 3))
