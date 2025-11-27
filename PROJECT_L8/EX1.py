import random



BASES = ['A', 'T', 'C', 'G']

MIN_HOST_LENGTH = 200
MAX_HOST_LENGTH = 400
NUM_TE_TO_INSERT = random.randint(3, 4)
TE_CORE_SEQUENCE = "GATTACA"  
IR_SEQUENCE = "TCGA"

TSD_LENGTH = 4
def generate_random_dna(length):
    
    return "".join(random.choice(BASES) for _ in range(length))

def get_complement(base):

    if base == 'A': return 'T'
    if base == 'T': return 'A'
    if base == 'C': return 'G'
    if base == 'G': return 'C'
    return base 

def get_reverse_complement(sequence):
    
    complement = "".join(get_complement(base) for base in sequence)
    return complement[::-1]


def create_simulated_dna_with_tes():
  
    host_length = random.randint(MIN_HOST_LENGTH, MAX_HOST_LENGTH)
    host_dna = list(generate_random_dna(host_length))
    te_positions = []
    
   
    ir2_sequence = get_reverse_complement(IR_SEQUENCE)
    
    
    TE_INSERT = IR_SEQUENCE + TE_CORE_SEQUENCE + ir2_sequence
    TE_INSERT_LENGTH = len(TE_INSERT)
    
    
    possible_indices = list(range(TSD_LENGTH, len(host_dna) - TSD_LENGTH))
    
  
    random.shuffle(possible_indices)

    print(f"--- 1. Generating DNA with {NUM_TE_TO_INSERT} TEs ---")
    print(f"Host DNA Length: {len(host_dna)}")
    
    inserted_count = 0
    
   
    for insertion_index in possible_indices:
        if inserted_count >= NUM_TE_TO_INSERT:
            break

        
        is_overlapping = False
        for start, end in te_positions:
            
            if (insertion_index > start - (TE_INSERT_LENGTH + TSD_LENGTH) and 
                insertion_index < end + (TE_INSERT_LENGTH + TSD_LENGTH)):
                is_overlapping = True
                break
        
        if is_overlapping:
            continue

        
        tsd_sequence = "".join(host_dna[insertion_index : insertion_index + TSD_LENGTH])

        
        full_insertion = tsd_sequence + TE_INSERT + tsd_sequence

       
        
        
        pre_insertion = host_dna[:insertion_index]
        
        post_insertion = host_dna[insertion_index + TSD_LENGTH:]

        
        host_dna = pre_insertion + list(full_insertion) + post_insertion

        
        te_start = len(pre_insertion) + len(tsd_sequence)
        te_end = te_start + TE_INSERT_LENGTH - 1
        te_positions.append((te_start, te_end))

        inserted_count += 1
        
        
        if inserted_count >= NUM_TE_TO_INSERT:
            break
        
    final_dna = "".join(host_dna)
    print(f"Final DNA Length: {len(final_dna)}")
    print(f"Actual TE Positions (0-indexed start, end): {te_positions}")
    print("-" * 50)
    
    return final_dna, te_positions, TE_INSERT, IR_SEQUENCE, ir2_sequence, TSD_LENGTH



def detect_transposable_elements(dna_sequence, te_insert_sequence, 
                                 ir_sequence_fwd, ir_sequence_rev, tsd_length):
    
    detected_positions = []
    
   
    te_length = len(te_insert_sequence)
    
    print("--- 2. Implementing Detection Algorithm ---")
    print(f"Searching for TE core sequence: {te_insert_sequence} (Length: {te_length})")
    print(f"TSD (Direct Repeat) Length: {tsd_length}")
    
   
    for i in range(len(dna_sequence) - te_length):
       
        potential_te = dna_sequence[i : i + te_length]
        
        if potential_te == te_insert_sequence:
           
            tsd1_start = i - tsd_length
            tsd1_end = i
            
            
            tsd2_start = i + te_length
            tsd2_end = i + te_length + tsd_length
            
            
            if tsd1_start < 0 or tsd2_end > len(dna_sequence):
                continue
                
            
            tsd1 = dna_sequence[tsd1_start : tsd1_end]
            tsd2 = dna_sequence[tsd2_start : tsd2_end]
            
            
            if tsd1 == tsd2:
                
                te_start = i
                te_end = i + te_length - 1
                detected_positions.append((te_start, te_end))
                
               

    print(f"Detected TE Positions (0-indexed start, end): {detected_positions}")
    print("-" * 50)
    return detected_positions



if __name__ == "__main__":
    
    final_dna, actual_te_positions, te_insert, ir1, ir2, tsd_len = create_simulated_dna_with_tes()

    
    detected_te_positions = detect_transposable_elements(
        final_dna, te_insert, ir1, ir2, tsd_len
    )
    
    
    print("### Comparison of Results ###")
    
    
    actual_te_positions.sort()
    detected_te_positions.sort()
    
    if actual_te_positions == detected_te_positions:
        print(" Detection Successful: It's a match!")
    else:
        print("Detection Failed:")
        print(f"   Actual Positions:   {actual_te_positions}")
        print(f"   Detected Positions: {detected_te_positions}")
    
    
    print("\nSnippet of DNA with the first detected TE highlighted:")
    if detected_te_positions:
        start, end = detected_te_positions[0]
        tsd_start = start - tsd_len
        tsd_end = end + tsd_len + 1 
        
  
        tsd1_seq = final_dna[tsd_start:start]
        te_seq = final_dna[start:end + 1]
        tsd2_seq = final_dna[end + 1:tsd_end]
        
        highlighted_snippet = (
            final_dna[tsd_start - 20:tsd_start] + 
            "" + tsd1_seq + "" +
            "---" + te_seq + "---" + 
            "" + tsd2_seq + "" + 
            final_dna[tsd_end:tsd_end + 20] 
        )
        
        print(f"Context... **TSD1---TE_CORE---TSD2** ...Context")
        print(highlighted_snippet)