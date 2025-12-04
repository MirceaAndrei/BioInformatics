import os
import math
import matplotlib.pyplot as plt

ENZYMES = {
    "EcoRI":  ("GAATTC", 1),
    "BamHI":  ("GGATCC", 1),
    "HindIII":("AAGCTT", 1),
    "TaqI":   ("TCGA",   1),
    "HaeIII": ("GGCC",   2),
}

def read_fasta(path):
    s = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            s.append(line.upper())
    return "".join(s)

def find_sites(seq, pat, cut):
    p = []
    i = 0
    while True:
        j = seq.find(pat, i)
        if j == -1:
            break
        p.append(j + cut)
        i = j + 1
    return p

def get_fragments(L, cuts):
    if not cuts:
        return [L]
    cuts = [0] + cuts + [L]
    fr = []
    for i in range(len(cuts)-1):
        fr.append(cuts[i+1] - cuts[i])
    return fr

def simulate_gel(data):
    allf = [v for x in data.values() for v in x]
    logs = [math.log10(v) for v in allf]
    mn, mx = min(logs), max(logs)
    sp = mx - mn if mx != mn else 1
    fig, ax = plt.subplots(figsize=(14,6))
    lane_s = 1.2
    lane_w = 0.5
    for i, (name, frs) in enumerate(data.items(), 1):
        xc = i * lane_s
        for L in frs:
            y = 1 - ((math.log10(L) - mn) / sp)
            ax.hlines(y, xc - lane_w/2, xc + lane_w/2, linewidth=4)
    ax.set_xlim(0, lane_s*(len(data)+1))
    ax.set_ylim(0,1)
    ax.invert_yaxis()
    ax.set_xticks([lane_s*i for i in range(1, len(data)+1)])
    ax.set_xticklabels(list(data.keys()), rotation=45)
    ax.set_yticklabels([])
    ax.set_ylabel("Normalized fragment migration")
    ax.set_title("Simulated Agarose Gel")
    plt.tight_layout()
    plt.show()

def main():
    fasta_files = [f for f in os.listdir() if f.lower().endswith(".fasta")]
    if not fasta_files:
        print("No FASTA file found.")
        return

    fasta = fasta_files[0]
    seq = read_fasta(fasta)
    L = len(seq)

    print("Sequence length:", L, "bp\n")

    results = {}

    for name, (pat, cut) in ENZYMES.items():
        sites = find_sites(seq, pat, cut)
        fr = get_fragments(L, sites)
        results[name] = fr
        print(name)
        print("Cuts:", [s+1 for s in sites])
        print("Fragments:", fr, "\n")

    simulate_gel(results)

main()
