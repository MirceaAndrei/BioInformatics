import os
from Bio import SeqIO
from Bio.Seq import Seq


MIN_IR_LENGTH = 4
MAX_IR_LENGTH = 6

MAX_TE_LENGTH = 5000 


def find_inverted_repeats(sequence_record):
   
    results = []
    seq = str(sequence_record.seq)
    genome_length = len(seq)
    
    print(f"-> Analyzing sequence: {sequence_record.id} (Length: {genome_length})")
    
  
    for i in range(genome_length):
        

        for ir_len in range(MIN_IR_LENGTH, MAX_IR_LENGTH + 1):
            
            
            if i + ir_len > genome_length:
                continue
                
            ir1_sequence = seq[i : i + ir_len]
            expected_ir2 = str(Seq(ir1_sequence).reverse_complement())
            start_search = i + ir_len
            end_search = min(genome_length, i + ir_len + MAX_TE_LENGTH)
            
           
            for j in range(start_search, end_search):
                
               
                if j + ir_len > genome_length:
                    continue
                    
                ir2_sequence = seq[j : j + ir_len]
                
               
                if ir2_sequence == expected_ir2:
                   
                    te_start = i + 1 
                    te_end = j + ir_len
                    te_length = te_end - te_start + 1
              
                    results.append({
                        'genome_id': sequence_record.id,
                        'position_start': te_start,
                        'position_end': te_end,
                        'length': te_length,
                        'ir_length': ir_len,
                        'ir1_seq': ir1_sequence
                    })
                    
    return results

def process_genomes(genome_files):
   
    all_results = []
    
    for file_path in genome_files:
        print(f"\n--- Processing File: {file_path} ---")
        
        
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                te_detections = find_inverted_repeats(record)
                all_results.extend(te_detections)
                
        except FileNotFoundError:
            print(f"ERROR: File not found at {file_path}. Skipping.")
        except Exception as e:
            print(f"An error occurred while processing {file_path}: {e}")
            
    return all_results

def display_results(results):
    
    if not results:
        print("\n--- RESULTS ---")
        print("No potential Transposable Elements detected based on criteria.")
        return

    print("\n--- DETECTED TRANSPOSEABLE ELEMENTS (Based on IR Structure) ---")
    print(f"{'Genome ID':<20} | {'Start':<8} | {'End':<8} | {'Length (bp)':<12} | {'IR Length':<10} | {'IR1 Sequence':<15}")
    print("-" * 80)
    for res in results:
        print(
            f"{res['genome_id'][:20]:<20} | {res['position_start']:<8} | {res['position_end']:<8} | {res['length']:<12} | {res['ir_length']:<10} | {res['ir1_seq']:<15}"
        )



if __name__ == "__main__":
    
  
    
  
    GENOME_PATHS = [
        "path/to/your/genome1.fasta",
        "path/to/your/genome2.fasta",
        "path/to/your/genome3.fasta",
    ]


    try:
        os.makedirs("test_data", exist_ok=True)
        TEST_FILE = "test_data/test_genome.fasta"
       
        test_seq = ">E_coli_K12\n" \
                   "ATATAC*TTAAGCTTTTTTTTTTCTTAAGATGATGGTCAGTTGTGTGTGTTGAC*AAAAA"
        with open(TEST_FILE, "w") as f:
            f.write(test_seq)
        GENOME_PATHS = [TEST_FILE] 
    except Exception as e:
        print(f"Could not create test file: {e}")

    detected_elements = process_genomes(GENOME_PATHS)
    display_results(detected_elements)