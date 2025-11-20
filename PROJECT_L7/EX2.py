import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import textwrap
import os
import sys

def read_multi_fasta(fasta_file):
    
    sequences = {}
    current_header = None
    current_sequence = []
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    if current_header and current_sequence:
                        sequences[current_header] = "".join(current_sequence).upper()
                    
                   
                    current_header = line[1:].split()[0]
                    current_sequence = []
                else:
                    current_sequence.append(line)
        
        
        if current_header and current_sequence:
            sequences[current_header] = "".join(current_sequence).upper()
            
    except FileNotFoundError:
        print(f"Error: The file '{fasta_file}' was not found.")
        return None
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None
        
    return sequences



def find_tandem_repeats(sequence, min_len=3, max_len=6):
   
    all_repeats = defaultdict(lambda: defaultdict(list))
    if not sequence:
        return all_repeats

    seq_len = len(sequence)
    
   
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
                    
                   
                    last_end_index = all_repeats[k][motif][-1][0] + k * all_repeats[k][motif][-1][1] - 1 if all_repeats[k][motif] else -1
                    
                    if start_index > last_end_index:
                        all_repeats[k][motif].append((start_index, repeat_count))
    
    return all_repeats

def analyze_genome_repeats(sequence):
  
    results = find_tandem_repeats(sequence)
    
    aggregated_data = []

    
    for k, motifs in results.items():
        for motif, locations in motifs.items():
            total_copies = sum(count for start, count in locations)
            aggregated_data.append((motif, total_copies))

    if not aggregated_data:
        return None, 0
    
    
    top_motif, max_copies = max(aggregated_data, key=lambda x: x[1])
    
    return top_motif, max_copies



def process_and_plot(fasta_file="sequences.fasta"):
    
    print(f" Reading and analyzing sequences from **{fasta_file}**...")
    genomes = read_multi_fasta(fasta_file)
    
    if not genomes:
        sys.exit(1)

    print(f"Successfully read {len(genomes)} genomes. Detecting top tandem repeats...")
    
    plot_data = []
    
   
    for header, sequence in genomes.items():
        top_motif, max_copies = analyze_genome_repeats(sequence)
        
       
        plot_id = textwrap.shorten(header, width=20, placeholder="...")

        plot_data.append({
            'Genome ID': plot_id,
            'Top Motif': top_motif if top_motif else 'None',
            'Total Copies': max_copies
        })
        
    df = pd.DataFrame(plot_data)
    
    if df.empty or df['Total Copies'].sum() == 0:
        print("\n✅ No significant tandem repeats (3-6 bases, min 2 copies) found in the genomes.")
        return

   
    df_sorted = df.sort_values(by='Total Copies', ascending=False).reset_index(drop=True)
    
    
    plt.figure(figsize=(14, 7))
    
   
    labels = [f"{row['Top Motif']} ({row['Total Copies']})" for index, row in df_sorted.iterrows()]
    
    plt.bar(df_sorted['Genome ID'], df_sorted['Total Copies'], color='#1f77b4')

   
    for i, (label, count) in enumerate(zip(labels, df_sorted['Total Copies'])):
      
        plt.text(i, count, label, ha='center', va='bottom', fontsize=9, rotation=0)

    plt.title('Most Frequent Tandem Repeats (3-6 bases, min 2 copies) in Genomes')
    plt.xlabel('Genome ID')
    plt.ylabel('Total Copies of Top Motif')
    plt.xticks(rotation=45, ha='right')
    plt.ylim(top=df_sorted['Total Copies'].max() * 1.2)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    plot_filename = 'most_frequent_repeats_chart.png'
    plt.savefig(plot_filename)
    print(f"\n🎉 Successfully generated plot: **{plot_filename}**")
    
    print("\nDataFrame of Top Repeats per Genome:")
    print(df_sorted[['Genome ID', 'Top Motif', 'Total Copies']].to_string(index=False))


if __name__ == "__main__":
    process_and_plot("sequences.fasta")