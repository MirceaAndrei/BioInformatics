
import math
import matplotlib.pyplot as plt

def read_fasta(filename):
    seqs = {}
    with open(filename, 'r') as f:
        name = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = ""
            else:
                seqs[name] += line.upper()
    return seqs

matrix = {
    'A': [0.3, 0.6, 0.1, 0.0, 0.0, 0.6, 0.7, 0.2, 0.1],
    'C': [0.2, 0.2, 0.1, 0.0, 0.0, 0.2, 0.1, 0.1, 0.2],
    'G': [0.1, 0.1, 0.7, 1.0, 0.0, 0.1, 0.1, 0.5, 0.1],
    'T': [0.4, 0.1, 0.1, 0.0, 1.0, 0.1, 0.1, 0.2, 0.6]
}

pseudo = 0.001
background = 0.25
pwm = {}
for nt, probs in matrix.items():
    pwm[nt] = [math.log((p + pseudo) / (background + 4 * pseudo)) for p in probs]

def get_scores(sequence):
    w = 9
    scores = []
    for i in range(len(sequence) - w + 1):
        window = sequence[i:i+w]
        s = sum(pwm[base][j] for j, base in enumerate(window) if base in pwm)
        scores.append(s)
    return scores

genomes = read_fasta("influenza.fasta")

for name, seq in genomes.items():
    scores = get_scores(seq)
    plt.figure(figsize=(12, 4))
    plt.plot(scores, color='blue', linewidth=0.5)
    plt.title(f"Motif Signal: {name}")
    plt.xlabel("Genome Position")
    plt.ylabel("Log-Likelihood Score")
    
    threshold = max(scores) * 0.8
    peaks = [i for i, s in enumerate(scores) if s > threshold]
    for p in peaks:
        plt.axvline(x=p, color='red', alpha=0.3, linestyle='--')
        
    plt.tight_layout()
    plt.show()

    print(f"Top signal for {name} found at position {scores.index(max(scores))}")