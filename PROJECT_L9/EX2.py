
import os
import math
import matplotlib.pyplot as plt

ENZYMES = {
    "EcoRI":  ("GAATTC", 1),
    "BamHI":  ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "TaqI":   ("TCGA", 1),
    "HaeIII": ("GGCC", 2),
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
    for i in range(len(cuts) - 1):
        fr.append(cuts[i + 1] - cuts[i])
    return fr

def digest_sequence(seq):
    L = len(seq)
    res = {}
    for name, (pat, cut) in ENZYMES.items():
        sites = find_sites(seq, pat, cut)
        fr = get_fragments(L, sites)
        res[name] = fr
    return res

def simulate_gel(lanes, title, save_path=None):
    allf = [v for x in lanes.values() for v in x]
    if not allf:
        print("No fragments to display for:", title)
        return
    logs = [math.log10(v) for v in allf]
    mn, mx = min(logs), max(logs)
    sp = mx - mn if mx != mn else 1
    width = max(10, 1.5 * len(lanes))
    fig, ax = plt.subplots(figsize=(width, 6))
    lane_s = 1.2
    lane_w = 0.5
    for i, (name, frs) in enumerate(lanes.items(), 1):
        xc = i * lane_s
        for L in frs:
            y = 1 - ((math.log10(L) - mn) / sp)
            ax.hlines(y, xc - lane_w / 2, xc + lane_w / 2, linewidth=4)
    ax.set_xlim(0, lane_s * (len(lanes) + 1))
    ax.set_ylim(0, 1)
    ax.invert_yaxis()
    ax.set_xticks([lane_s * i for i in range(1, len(lanes) + 1)])
    ax.set_xticklabels(list(lanes.keys()), rotation=45)
    ax.set_yticklabels([])
    ax.set_ylabel("Normalized fragment migration")
    ax.set_title(title)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.show()
    plt.close(fig)

def main():
    fasta_files = [f for f in os.listdir() if f.lower().endswith(".fasta")]
    fasta_files = [f for f in fasta_files if "covid" not in f.lower() and "ebola" not in f.lower()]
    fasta_files.sort()
    if len(fasta_files) == 0:
        print("No influenza FASTA files found.")
        return
    if len(fasta_files) > 10:
        fasta_files = fasta_files[:10]

    print("Genomes used:")
    for f in fasta_files:
        print(" -", f)

    genome_digests = {}

    for fasta in fasta_files:
        seq = read_fasta(fasta)
        print("\n", fasta, "length:", len(seq), "bp")
        digest = digest_sequence(seq)
        genome_digests[fasta] = digest
        for enz, fr in digest.items():
            print(enz, "fragments:", fr)
        simulate_gel(digest, title=f"{fasta} – Restriction Digest Gel", save_path=f"{fasta}_gel.png")

    common_per_enzyme = {}
    for enz in ENZYMES.keys():
        sets = [set(genome_digests[g][enz]) for g in fasta_files]
        if sets:
            inter = sets[0]
            for s in sets[1:]:
                inter = inter.intersection(s)
            common_per_enzyme[enz] = inter
        else:
            common_per_enzyme[enz] = set()

    diffs_per_genome_per_enzyme = {}
    for g in fasta_files:
        diffs_per_genome_per_enzyme[g] = {}
        for enz in ENZYMES.keys():
            fr = genome_digests[g][enz]
            common = common_per_enzyme[enz]
            diffs = [x for x in fr if x not in common]
            diffs_per_genome_per_enzyme[g][enz] = diffs
            print(g, enz, "unique fragments:", diffs)
        simulate_gel(diffs_per_genome_per_enzyme[g],
                     title=f"{g} – Unique Fragments (Common Bands Removed)",
                     save_path=f"{g}_unique_gel.png")

    merged_lanes = {}
    for g in fasta_files:
        all_unique = []
        for enz in ENZYMES.keys():
            all_unique.extend(diffs_per_genome_per_enzyme[g][enz])
        merged_lanes[g] = all_unique

    simulate_gel(merged_lanes,
                 title="Combined Gel – Unique Fragments Across All Genomes",
                 save_path="combined_unique_fragments_gel.png")

main()
