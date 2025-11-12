
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

# 2. Transcription Function (DNA -> RNA)
def transcribe(dna_sequence):
    """
    Converts a DNA coding strand sequence into an mRNA sequence.
    Replaces all 'T' (Thymine) with 'U' (Uracil).
    """
    # Use .upper() to make the function case-insensitive
    return dna_sequence.upper().replace('T', 'U')

# 3. Translation Function (RNA -> Amino Acid Sequence)
def translate(rna_sequence):
    """
    Converts an mRNA sequence into an amino acid (protein) sequence.
    Translation starts at the first 'AUG' (Met) codon and stops
    when a 'STOP' codon is encountered.
    """
    protein = []
    is_translating = False

    # Find the first start codon 'AUG'
    start_index = rna_sequence.find('AUG')
    
    if start_index == -1:
        # If no start codon is found, translation cannot begin
        return "No start codon (AUG) found in sequence."

    # Iterate over the RNA sequence in steps of 3 (one codon at a time)
    # Start iterating from the 'AUG' codon
    for i in range(start_index, len(rna_sequence), 3):
        # Get the 3-letter codon
        codon = rna_sequence[i:i+3]

        # Ensure we have a full codon
        if len(codon) < 3:
            break  # Reached the end of the sequence

        # Look up the codon in the genetic code table
        # .get(codon, 'X') returns 'X' (unknown) if the codon is not found
        amino_acid = GENETIC_CODE.get(codon, 'X')

        # Check for stop codon
        if amino_acid == 'STOP':
            break  # Stop translation
        
        # Add the amino acid to our protein sequence
        protein.append(amino_acid)

    # Join the list of amino acids into a single string, separated by dashes
    return "-".join(protein)

# --- Main Application ---
if __name__ == "__main__":
    # Example DNA coding sequence (from ATG to TAA)
    # This is a simplified gene sequence for demonstration
    dna_coding_strand = "GGCATGTACCCGGATTGTCGTCAAGCCGCTAGCATAA"

    print(f"Original DNA Sequence: {dna_coding_strand}")
    print("-" * 30)

    # Step 1: Transcribe DNA to RNA
    rna_sequence = transcribe(dna_coding_strand)
    print(f"Transcribed RNA Sequence: {rna_sequence}")
    print("-" * 30)

    # Step 2: Translate RNA to Protein
    protein_sequence = translate(rna_sequence)
    print(f"Translated Protein Sequence: {protein_sequence}")