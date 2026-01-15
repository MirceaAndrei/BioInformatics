
import numpy as np
import matplotlib.pyplot as plt
import io

WINDOW_LENGTH = 500
COV_COLOR = 'darkblue'
FLU_COLOR = 'darkred'

def read_fasta(filename):
    sequences = []
    current_header = None
    current_sequence = []
    
    try:
        with open(filename, 'r') as f:
            data = f.read()
    except FileNotFoundError:
        print(f"Warning: File '{filename}' not found. Returning empty list of sequences.")
        return sequences
            
    for line in data.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_header and current_sequence:
                sequences.append({'header': current_header, 'sequence': "".join(current_sequence)})
            current_header = line[1:].split('|')[0]
            current_sequence = []
        else:
            current_sequence.append(line.upper().replace('N', ''))

    if current_header and current_sequence:
        sequences.append({'header': current_header, 'sequence': "".join(current_sequence)})
        
    return sequences

def calculate_cg_percent(sequence):
    N = len(sequence)
    if N == 0:
        return 0.0
    C_count = sequence.count('C')
    G_count = sequence.count('G')
    return ((C_count + G_count) / N) * 100.0

def run_sliding_window_analysis(sequence, window):
    cg_values = []
    window_centers = []
    
    for i in range(0, len(sequence) - window + 1, window // 5):
        window_sequence = sequence[i:i + window]
        cg_values.append(calculate_cg_percent(window_sequence))
        window_centers.append(i + window / 2)
        
    return np.array(window_centers), np.array(cg_values)

def calculate_center_of_weight(positions, values):
    numerator = np.sum(positions * values)
    denominator = np.sum(values)
    if denominator == 0:
        return 0
    return numerator / denominator

def plot_objective_digital_straints(all_data):
    plt.figure(figsize=(15, 8))
    
    for data in all_data:
        label_base = data['header'].split(' ')[0]
        color = COV_COLOR if 'SARS-CoV-2' in data['header'] or 'Severe acute' in data['header'] else FLU_COLOR
        plt.plot(data['positions'], data['cg_values'], label=data['short_label'], color=color, alpha=0.7, linewidth=1.5)

    plt.title('Objective Digital Straint (C+G % in 500bp Window) for Viral Genomes')
    plt.xlabel('Genomic Position (bp)')
    plt.ylabel('C+G %')
    plt.legend(ncol=4, fontsize='small', loc='upper right')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()

def plot_centers_of_weight(all_data):
    plt.figure(figsize=(15, 10))
    
    cow_values = []
    avg_cg_values = []
    labels = []
    colors = []
    
    for data in all_data:
        cow_values.append(data['cow'])
        avg_cg_values.append(data['avg_cg'])
        labels.append(data['short_label'])
        colors.append(COV_COLOR if 'SARS-CoV-2' in data['header'] or 'Severe acute' in data['header'] else FLU_COLOR)

    scatter = plt.scatter(cow_values, avg_cg_values, s=150, c=colors, alpha=0.8, edgecolors='black', linewidth=1, zorder=3)
    
    for i, label in enumerate(labels):
        plt.annotate(label, (cow_values[i], avg_cg_values[i]), 
                     textcoords="offset points", xytext=(5, 5), 
                     ha='left', fontsize=8)

    plt.title('Center of Weight Distribution for Objective Digital Straints')
    plt.xlabel('Center of Weight Position (bp)')
    plt.ylabel('Average C+G % (Full Sequence)')
    plt.grid(True, linestyle='--', alpha=0.5)
    
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label='COVID-19', markerfacecolor=COV_COLOR, markersize=10),
                       plt.Line2D([0], [0], marker='o', color='w', label='Influenza', markerfacecolor=FLU_COLOR, markersize=10)]
    plt.legend(handles=legend_elements, loc='lower right')
    plt.show()

all_genome_data = []

covid_sequences = read_fasta('covid.fasta')
flu_sequences = read_fasta('sequences.fasta')

all_sequences = covid_sequences + flu_sequences

for seq_data in all_sequences:
    header = seq_data['header']
    sequence = seq_data['sequence']
    
    species = 'COVID-19' if 'SARS-CoV-2' in header or 'Severe acute' in header else 'Flu'
    short_label = f"{species}: {header.split(' ')[-1]}"
    
    if len(sequence) < WINDOW_LENGTH:
        continue

    positions, cg_values = run_sliding_window_analysis(sequence, WINDOW_LENGTH)
    cow = calculate_center_of_weight(positions, cg_values)
    avg_cg = calculate_cg_percent(sequence)
    
    all_genome_data.append({
        'header': header,
        'short_label': short_label,
        'positions': positions,
        'cg_values': cg_values,
        'cow': cow,
        'avg_cg': avg_cg
    })

plot_objective_digital_straints(all_genome_data)
plot_centers_of_weight(all_genome_data)