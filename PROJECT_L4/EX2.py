import sys
import os
from collections import Counter
import matplotlib.pyplot as plt
# Note: urllib and json imports have been removed.

# 1. The Genetic Code Table (from your image)
GENETIC_CODE = {
    # U-block
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
    # C-block
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    # A-block
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    # G-block
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

# 2. Pre-loaded Food Information Database
AMINO_ACID_INFO = {
    'Leu': "Leucine: An essential amino acid. Found in high-protein foods (meat, dairy, soy). Some grains like corn are lower in leucine.",
    'Ser': "Serine: A non-essential amino acid, meaning the body can produce it. It's abundant in many foods.",
    'Thr': "Threonine: An essential amino acid. Often found in high-protein sources. Grains like wheat and rice can be lower in threonine.",
    'Arg': "Arginine: A conditionally-essential amino acid. Abundant in nuts, seeds, and meats. Dairy is a source, but generally lower than nuts.",
    'Lys': "Lysine: An essential amino acid. Famously low in grains (e.g., corn, wheat, rice). Abundant in legumes (beans, lentils), meat, and dairy.",
    'Glu': "Glutamic acid (Glutamate): A non-essential amino acid. Extremely common and found in almost all protein-containing foods (e.g., tomatoes, cheese, mushrooms).",
    'Ala': "Alanine: A non-essential amino acid. Abundant in meat, poultry, fish, and dairy products.",
    'Asn': "Asparagine: A non-essential amino acid. Found in dairy, beef, poultry, and potatoes.",
    'Gln': "Glutamine: A non-essential amino acid. Abundant in both animal and plant proteins, including cabbage and beets.",
    'Val': "Valine: An essential amino acid. Found in high-protein foods, soy, and peanuts.",
    'Phe': "Phenylalanine: An essential amino acid. Found in most protein-rich foods, especially soy, nuts, and seeds."
}

# General list of foods that are inherently low in all proteins
LOW_PROTEIN_FOODS = {
    "Fats and Oils (Typically protein-free)": [
        "Olive oil, coconut oil, vegetable oils",
        "Butter, margarine, mayonnaise"
    ],
    "Sugars and Simple Starches (Primarily carbohydrates)": [
        "Table sugar, corn-starch, honey, maple syrup",
        "Sorbets, hard candies, jams, jellies"
    ],
    "Certain Fruits and Vegetables (Very low in protein)": [
        "Apples, grapes, berries",
        "Cucumber, celery, lettuce, bell peppers, carrots"
    ],
    "Beverages": [
        "Water, coffee, tea (without milk/cream)",
        "Most fruit juices (apple, grape)"
    ]
}


def parse_fasta(filename):
    """
    Parses a FASTA file and returns the complete genome sequence as a
    single, uppercase string.
    """
    if not os.path.exists(filename):
        print(f"Error: File not found at '{filename}'", file=sys.stderr)
        print("Please make sure the file is in the same folder as the script.", file=sys.stderr)
        sys.exit(1)
        
    print(f"Parsing {filename}...")
    sequence = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            sequence.append(line.strip())
            
    return "".join(sequence).upper()

def transcribe(dna_sequence):
    """
    Converts a DNA coding strand sequence into an mRNA sequence.
    """
    return dna_sequence.replace('T', 'U')

def get_codon_frequencies(rna_sequence):
    """
    Counts the frequency of each 3-letter codon in an RNA sequence.
    """
    codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence) - 2, 3)]
    return Counter(codons)

def get_amino_acid_frequencies(codon_counts, genetic_code):
    """
    Calculates the total frequency of each amino acid based on codon counts.
    """
    amino_acid_counts = Counter()
    for codon, count in codon_counts.items():
        amino_acid = genetic_code.get(codon)
        if amino_acid and amino_acid != 'STOP':
            amino_acid_counts[amino_acid] += count
    return amino_acid_counts

def plot_top_codons(codon_counts, title, top_n=10):
    """
    Creates a bar chart for the top N most frequent codons.
    """
    top_codons = codon_counts.most_common(top_n)
    codons, counts = zip(*top_codons)
    
    plt.figure(figsize=(12, 7))
    plt.bar(codons, counts, color='darkgreen')
    plt.title(title, fontsize=16)
    plt.xlabel("Codon", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.tight_layout()

def print_top_amino_acids(amino_acid_counts, title, top_n=3):
    """
    Prints the top N most frequent amino acids to the console.
    """
    print(f"\n>>> TOP 3 AMINO ACIDS ({title})")
    for i, (aa, count) in enumerate(amino_acid_counts.most_common(top_n), 1):
        print(f"{i}. {aa}: {count} occurrences")

def find_and_print_food_recommendations(amino_acid_set):
    """
    Prints food recommendations from the pre-loaded database based on
    the set of top amino acids.
    """
    print("\n" + "="*50)
    print("STATIC FOOD RECOMMENDATION (Pre-loaded)")
    print("="*50)
    
    amino_acid_list_str = ", ".join(sorted(list(amino_acid_set)))
    print(f"The top amino acids from both genomes are: {amino_acid_list_str}.")
    print("Here is a breakdown of those amino acids and a general list of")
    print("foods that are low in protein (and thus low in these amino acids).\n")
    
    print("--- Notes on Top Amino Acids ---")
    # Loop through the set of top amino acids and print info
    for aa in sorted(list(amino_acid_set)):
        # Find the info from our database, or provide a default message
        info = AMINO_ACID_INFO.get(aa, f"{aa}: No specific info pre-loaded.")
        print(f"* {info}")
    
    print("\n--- General Low-Protein Food Categories ---")
    print("To avoid all these amino acids, you would seek low-protein foods:")
    
    # Loop through the low-protein food database and print it
    for category, items in LOW_PROTEIN_FOODS.items():
        print(f"\n{category}:")
        for item in items:
            print(f"  - {item}")
    
    print("="*50)


def main():
    """
    Main function to run the complete analysis.
    """
    COVID_FASTA = "covid.fasta"
    INFLUENZA_FASTA = "influenza.fasta"

    # --- COVID-19 Analysis ---
    covid_dna = parse_fasta(COVID_FASTA)
    covid_rna = transcribe(covid_dna)
    covid_codon_counts = get_codon_frequencies(covid_rna)
    covid_aa_counts = get_amino_acid_frequencies(covid_codon_counts, GENETIC_CODE)
    plot_top_codons(covid_codon_counts, "Top 10 Most Frequent Codons: COVID-19")

    # --- Influenza Analysis ---
    flu_dna = parse_fasta(INFLUENZA_FASTA)
    flu_rna = transcribe(flu_dna)
    flu_codon_counts = get_codon_frequencies(flu_rna)
    flu_aa_counts = get_amino_acid_frequencies(flu_codon_counts, GENETIC_CODE)
    plot_top_codons(flu_codon_counts, "Top 10 Most Frequent Codons: Influenza")

    # --- Console Output ---
    print("\n" + "="*50)
    print("GENOME ANALYSIS CONSOLE OUTPUT")
    print("="*50)

    # Find common codons
    covid_top_10 = set([codon for codon, count in covid_codon_counts.most_common(10)])
    flu_top_10 = set([codon for codon, count in flu_codon_counts.most_common(10)])
    common_codons = covid_top_10.intersection(flu_top_10)
    
    print(f"\n--- Common Codons in Top 10 of Both Genomes ---")
    if common_codons:
        print(f"{'Codon':<6} | {'COVID-19 Count':<15} | {'Influenza Count':<15}")
        print("-" * 40)
        common_codon_stats = []
        for codon in common_codons:
            common_codon_stats.append((codon, covid_codon_counts[codon], flu_codon_counts[codon]))
        common_codon_stats.sort(key=lambda x: x[1], reverse=True)
        for codon, covid_count, flu_count in common_codon_stats:
            print(f"{codon:<6} | {covid_count:<15} | {flu_count:<15}")
    else:
        print("No codons were found in the top 10 of both genomes.")

    # --- Print Top 3 Amino Acids (as requested) ---
    print_top_amino_acids(covid_aa_counts, "COVID-19")
    print_top_amino_acids(flu_aa_counts, "Influenza")

    # --- Dynamic Food Recommendation ---
    covid_top_3_aa = set([aa for aa, count in covid_aa_counts.most_common(3)])
    flu_top_3_aa = set([aa for aa, count in flu_aa_counts.most_common(3)])
    all_top_aa = covid_top_3_aa.union(flu_top_3_aa)
    
    find_and_print_food_recommendations(all_top_aa)

    # --- Show Plots ---
    print("\nGenerating charts...")
    print("Please close the chart windows to exit the program.")
    plt.show()

if __name__ == "__main__":
    main()
