import sys
import os
from collections import defaultdict

def extract_sequence(fasta_file):

    try:
        with open(fasta_file, 'r') as f:
            lines = [line.strip() for line in f if not line.startswith('>')]
        
       
        if not lines:
            print(f"Error: No sequence found in {fasta_file}. The file might be empty or improperly formatted.")
            return None
            
        
        sequence = "".join(lines).upper()
        
       
        if any(base not in 'ATCG' for base in sequence):
            print("Warning: The sequence contains non-standard DNA bases. Proceeding, but results might be unexpected.")

        return sequence

    except FileNotFoundError:
        print(f"Error: The file '{fasta_file}' was not found.")
        return None
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None

def find_tandem_repeats(sequence, min_len=3, max_len=6):
   

    all_repeats = defaultdict(lambda: defaultdict(list))
    
    if not sequence:
        return all_repeats

    seq_len = len(sequence)
    print(f"\nSequence length: {seq_len} nucleotides.")
    
   
    for k in range(min_len, max_len + 1):
        
        for i in range(seq_len - k + 1):
           
            motif = sequence[i:i + k]
            
           
            next_start_index = i + k
            
    
            if next_start_index + k <= seq_len:
                next_motif = sequence[next_start_index:next_start_index + k]
                
               
                if motif == next_motif:
                    
                    repeat_count = 2 
                    current_end = next_start_index + k
                    
                    
                    while current_end + k <= seq_len and sequence[current_end:current_end + k] == motif:
                        repeat_count += 1
                        current_end += k
                    
                   
                    start_index = i
                    
                 
                    if not all_repeats[k][motif] or start_index > all_repeats[k][motif][-1][0]:
                        
                        
                        all_repeats[k][motif].append((start_index, repeat_count))
    
    return all_repeats

def print_results(results):
    
    has_repeats = False
    for k, motifs in sorted(results.items()):
        if motifs:
            has_repeats = True
            print(f"\n--- 📏 Motifs of Length {k} ---")
            
            
            for motif, locations in sorted(motifs.items()):
                print(f"**Motif: {motif}**")
                for start, count in locations:
                    end = start + (k * count) - 1
                    print(f"  - Location: {start} to {end} (0-indexed). Repeats: {count} copies. Total length: {k*count}b")
                print("-" * 20) 

    if not has_repeats:
        print("\n✅ No tandem repeats of lengths 3 to 6 (min 2 copies) were found in the sequence.")
    else:
        print("\n--- Detection Complete ---")



if __name__ == "__main__":
    FASTA_FILE_NAME = "covid.fasta"
    
    
    dna_sequence = extract_sequence(FASTA_FILE_NAME)
    
    if dna_sequence:
       
        repeat_data = find_tandem_repeats(dna_sequence)
        
        
        print_results(repeat_data)