import math

sequences = [
    "GAGGTAAAC", "TCCGTAGGT", "CAGGTTGGA", "ACAGT CAGT", 
    "TAGGTCATT", "TAGGTACTG", "ATGGT AACT", "CAGGTATAC", 
    "TGTGTGAGT", "AAGGTAAGT"
]

sequences = [s.replace(" ", "") for s in sequences]
n_seq = len(sequences)
seq_len = len(sequences[0])
bases = ['A', 'C', 'G', 'T']

count_matrix = {base: [0] * seq_len for base in bases}
for seq in sequences:
    for i, base in enumerate(seq):
        count_matrix[base][i] += 1

pseudo = 0.1
relative_freq_matrix = {base: [] for base in bases}
for base in bases:
    for i in range(seq_len):
        freq = (count_matrix[base][i] + pseudo) / (n_seq + 4 * pseudo)
        relative_freq_matrix[base].append(freq)

null_model = 0.25
log_likelihood_matrix = {base: [] for base in bases}
for base in bases:
    for freq in relative_freq_matrix[base]:
        log_likelihood_matrix[base].append(math.log(freq / null_model))

print("--- Count Matrix ---")
for b in bases: print(f"{b}: {count_matrix[b]}")

print("\n--- Relative Frequencies Matrix ---")
for b in bases: print(f"{b}: {[round(f, 3) for f in relative_freq_matrix[b]]}")

print("\n--- Log-likelihoods Matrix ---")
for b in bases: print(f"{b}: {[round(l, 3) for l in log_likelihood_matrix[b]]}")