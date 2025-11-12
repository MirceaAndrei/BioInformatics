import random
from Bio import Entrez, SeqIO
from io import StringIO
from collections import defaultdict, Counter
import matplotlib.pyplot as plt

Entrez.email = "your.email@example.com"  # <-- replace with your real email

# 1. Fetch a real DNA sequence from NCBI
def fetch_dna_sequence(min_len=1000, max_len=3000):
    handle = Entrez.esearch(db="nucleotide", term="human[organism] AND 1000:3000[SLEN]", retmax=5)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        raise ValueError("No sequences found.")

    for seq_id in record["IdList"]:
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
        seq_text = handle.read()
        handle.close()
        seq_record = SeqIO.read(StringIO(seq_text), "fasta")
        sequence = str(seq_record.seq).upper()
        seq_len = len(sequence)
        if min_len <= seq_len <= max_len:
            return sequence, seq_id, seq_record.description

    raise ValueError("No suitable sequence in range.")


# 2. Create random samples (reads)
def generate_reads(dna, num_reads=200, min_len=100, max_len=150):
    reads = []
    for _ in range(num_reads):
        read_length = random.randint(min_len, max_len)
        start = random.randint(0, len(dna) - read_length)
        reads.append((start, start + read_length))
    return reads


# 3. Visualize reads over original DNA sequence
def plot_reads(dna_length, reads, coverage):
    plt.figure(figsize=(12, 6))

    # Plot the original sequence line
    plt.hlines(y=0, xmin=0, xmax=dna_length, color='black', linewidth=2, label='Original DNA')

    # Plot each read as a horizontal colored segment
    for i, (start, end) in enumerate(reads):
        plt.hlines(y=i + 1, xmin=start, xmax=end,
                   color=plt.cm.tab20(i % 20), linewidth=3)

    plt.xlabel("DNA base position")
    plt.ylabel("Read index")
    plt.title("DNA Sequence Coverage Visualization")

    plt.tight_layout()
    plt.show()

    # Plot coverage depth
    plt.figure(figsize=(12, 3))
    plt.plot(range(dna_length), coverage, color='steelblue')
    plt.xlabel("DNA base position")
    plt.ylabel("Coverage depth")
    plt.title("Coverage across DNA sequence")
    plt.tight_layout()
    plt.show()


# 4. Reconstruct sequence by majority voting
def reconstruct_sequence(sequence, reads):
    pos_to_bases = defaultdict(list)
    for start, end in reads:
        for i in range(start, end):
            if i < len(sequence):
                pos_to_bases[i].append(sequence[i])  # Perfect sampling (no errors)

    reconstructed = []
    coverage = []
    for pos in range(len(sequence)):
        bases = pos_to_bases[pos]
        if not bases:
            reconstructed.append("N")
            coverage.append(0)
        else:
            count = Counter(bases)
            majority = count.most_common(1)[0][0]
            reconstructed.append(majority)
            coverage.append(len(bases))

    reconstructed_seq = "".join(reconstructed)
    errors = sum(a != b for a, b in zip(sequence, reconstructed_seq))
    accuracy = (len(sequence) - errors) / len(sequence) * 100
    avg_cov = sum(coverage) / len(coverage)

    return reconstructed_seq, accuracy, errors, coverage, avg_cov


# 5. Main program
if __name__ == "__main__":
    dna, seq_id, description = fetch_dna_sequence()

    print(f"Sequence ID: {seq_id}")
    print(f"Description: {description}")
    print(f"Length: {len(dna)} bases\n")

    reads = generate_reads(dna, 200, 100, 150)
    reconstructed, accuracy, errors, coverage, avg_cov = reconstruct_sequence(dna, reads)

    print(f"Reconstructed length: {len(reconstructed)}")
    print(f"Accuracy: {accuracy:.2f}%")
    print(f"Errors: {errors}")
    print(f"Average coverage: {avg_cov:.1f}x")

    print("\nMain problem: overlaps create conflicts.")
    print("High coverage gives >99% accuracy.")
    print("Real assembly needs graphs for repeats.")

    # Visualization (like your photo)
    plot_reads(len(dna), reads, coverage)
