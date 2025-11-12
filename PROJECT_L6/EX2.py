import random
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re
import os


def read_multi_fasta(filepath):
    if not os.path.exists(filepath):
        print(f"--- ERROR ---")
        print(f"File not found: {filepath}")
        print("Please make sure your FASTA file is in the same directory")
        print("and is named 'influenza_genomes.fasta'.")
        print("---------------")
        return None

    print(f"Reading sequences from {filepath}...")
    sequences = {}
    current_seq = ""
    current_header = ""

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                if current_header and current_seq:
                    sequences[current_header] = current_seq
                
                current_header = line[1:].split(' ')[0] 
                current_seq = ""
            else:
                current_seq += line
        
        if current_header and current_seq:
            sequences[current_header] = current_seq
            
    print(f"Successfully read {len(sequences)} sequences.")
    return sequences

def restriction_digest(sequence, enzyme_site='GAATTC'):
    cut_point = 1 
    
    matches = [m.start() for m in re.finditer(enzyme_site, sequence, re.IGNORECASE)]
    
    if not matches:
        return [len(sequence)]
    
    cut_indices = [m + cut_point for m in matches]
    
    fragment_lengths = []
    last_cut = 0
    
    for cut in cut_indices:
        fragment_lengths.append(len(sequence[last_cut:cut]))
        last_cut = cut
    
    fragment_lengths.append(len(sequence[last_cut:]))
    
    return fragment_lengths


def get_normalized_y_pos(length, log_min_all, log_max_all):
    if log_max_all == log_min_all: 
        return 0.5
    
    if length < 10: length = 10 
    
    log_len = math.log10(length)
    
    norm_pos = (log_len - log_min_all) / (log_max_all - log_min_all)
    norm_pos = max(0, min(1, norm_pos))
    return norm_pos


def plot_combined_gel(digests_dict, ladder_sizes):
    print("Generating combined gel plot...")
    
    num_lanes = len(digests_dict) + 1 
    fig, ax = plt.subplots(figsize=(max(12, num_lanes * 1.2), 10))
    
    LANE_WIDTH = 0.8
    BAND_HEIGHT = 0.03
    LANE_SPACING = 1.0

    all_lengths = list(ladder_sizes)
    for fragments in digests_dict.values():
        all_lengths.extend(fragments)
    min_len = min(all_lengths)
    max_len = max(all_lengths)
    if min_len < 10: min_len = 10
    
    log_min_all = math.log10(min_len)
    log_max_all = math.log10(max_len)

    ax.set_facecolor('black')
    ax.set_xlim(0, num_lanes * LANE_SPACING + 0.5)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    for i in range(num_lanes):
        x_center = (i * LANE_SPACING) + LANE_SPACING/2
        well = patches.Rectangle((x_center - LANE_WIDTH/2, 1.05),
                                 LANE_WIDTH, 0.05,
                                 facecolor='gray', edgecolor='white', lw=1.5, clip_on=False)
        ax.add_patch(well)
        
    x_ladder_center = LANE_SPACING / 2
    ax.text(x_ladder_center, 0.98, 'Ladder', ha='center', va='top', color='gray', fontsize=9)
    for size in ladder_sizes:
        y_pos = get_normalized_y_pos(size, log_min_all, log_max_all)
        ax.barh(y=y_pos, width=LANE_WIDTH, height=BAND_HEIGHT,
                left=x_ladder_center - LANE_WIDTH/2, color='white',
                edgecolor='darkgray', lw=0.5)
        ax.text(x_ladder_center - LANE_WIDTH/2 - 0.1, y_pos,
                f'{size} bp -', ha='right', va='center', color='black', fontsize=9,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, boxstyle='round,pad=0.2'))

    for i, (header, fragment_list) in enumerate(digests_dict.items()):
        x_sample_center = (i + 1) * LANE_SPACING + LANE_SPACING / 2
        short_header = header.replace('_', ' ').replace('Segment', 'Seg')
        ax.text(x_sample_center, 0.98, short_header, 
                ha='center', va='top', color='gray', fontsize=9)
        
        for length in fragment_list:
            y_pos = get_normalized_y_pos(length, log_min_all, log_max_all)
            ax.barh(y=y_pos, width=LANE_WIDTH, height=BAND_HEIGHT,
                    left=x_sample_center - LANE_WIDTH/2, color='white',
                    edgecolor='darkgray', lw=0.5)

    ax.set_title('Combined Gel Simulation: EcoRI Digest of 10 Influenza Segments', fontsize=16, color='black', pad=30)
    ax.text(num_lanes * LANE_SPACING / 2, 1.15, 'Wells (Negative Electrode)',
            ha='center', va='bottom', fontsize=10, color='gray')
    ax.text(num_lanes * LANE_SPACING / 2, -0.05, 'Positive Electrode',
            ha='center', va='top', fontsize=10, color='gray')

    plt.tight_layout()
    plt.savefig('combined_gel_simulation.png', facecolor='white')
    print("Combined gel plot saved as combined_gel_simulation.png")


def plot_separate_gels(digests_dict, ladder_sizes):
    print("Generating separate-lane gel plot...")
    
    num_segments = len(digests_dict)
    fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(10, 20), sharey=True)
    fig.suptitle('Separate Gels: EcoRI Digest of Each Segment', fontsize=20, y=1.03)
    
    ax_list = axes.flatten()

    all_lengths = list(ladder_sizes)
    for fragments in digests_dict.values():
        all_lengths.extend(fragments)
    min_len = min(all_lengths)
    max_len = max(all_lengths)
    if min_len < 10: min_len = 10
    
    log_min_all = math.log10(min_len)
    log_max_all = math.log10(max_len)

    LANE_WIDTH = 0.8
    BAND_HEIGHT = 0.03
    LANE_SPACING = 1.0
    NUM_LANES_PER_PLOT = 2
    
    segment_items = list(digests_dict.items())

    for i in range(len(ax_list)):
        ax = ax_list[i]
        
        if i >= len(segment_items):
            ax.axis('off') 
            continue
            
        header, fragment_list = segment_items[i]
        
        ax.set_facecolor('black')
        ax.set_xlim(0, NUM_LANES_PER_PLOT * LANE_SPACING + 0.5)
        ax.set_ylim(-0.1, 1.1)
        ax.set_xticks([])
        ax.set_yticks([]) 
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        
        ax.set_title(header.replace('_', ' '), fontsize=12, pad=15)
        
        for j in range(NUM_LANES_PER_PLOT):
            x_center = (j * LANE_SPACING) + LANE_SPACING/2
            well = patches.Rectangle((x_center - LANE_WIDTH/2, 1.05),
                                     LANE_WIDTH, 0.05,
                                     facecolor='gray', edgecolor='white', lw=1.5, clip_on=False)
            ax.add_patch(well)
        
        x_ladder_center = LANE_SPACING / 2
        ax.text(x_ladder_center, 0.98, 'Ladder', ha='center', va='top', color='gray', fontsize=9)
        for size in ladder_sizes:
            y_pos = get_normalized_y_pos(size, log_min_all, log_max_all)
            ax.barh(y=y_pos, width=LANE_WIDTH, height=BAND_HEIGHT,
                    left=x_ladder_center - LANE_WIDTH/2, color='white',
                    edgecolor='darkgray', lw=0.5)
            if i == 0:
                ax.text(x_ladder_center - LANE_WIDTH/2 - 0.1, y_pos,
                        f'{size} bp -', ha='right', va='center', color='black', fontsize=9,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, boxstyle='round,pad=0.2'))

        x_sample_center = 1 * LANE_SPACING + LANE_SPACING / 2
        ax.text(x_sample_center, 0.98, 'Sample', ha='center', va='top', color='gray', fontsize=9)
        for length in fragment_list:
            y_pos = get_normalized_y_pos(length, log_min_all, log_max_all)
            ax.barh(y=y_pos, width=LANE_WIDTH, height=BAND_HEIGHT,
                    left=x_sample_center - LANE_WIDTH/2, color='white',
                    edgecolor='darkgray', lw=0.5)
                    
    fig.text(0.5, 0.01, 'Positive Electrode', ha='center', va='bottom', fontsize=12, color='gray')
    fig.text(0.5, 0.98, 'Wells (Negative Electrode)', ha='center', va='top', fontsize=12, color='gray')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.97]) 
    plt.savefig('separate_gels_simulation.png', facecolor='white')
    print("Separate-lane gel plot saved as separate_gels_simulation.png")


if __name__ == "__main__":
    
    FASTA_FILE = 'influenza_genomes.fasta'
    ENZYME_SITE = 'GAATTC' 
    LADDER_SIZES = [3000, 2500, 2000, 1500, 1200, 1000, 800, 600, 500, 400, 300, 200, 100]

    
    sequences_dict = read_multi_fasta(FASTA_FILE)
    
    if sequences_dict:
        all_digests = {}
        print("\n--- 🧬 Performing In-Silico EcoRI Digest ---")
        for header, seq in sequences_dict.items():
            fragment_lengths = restriction_digest(seq, ENZYME_SITE)
            all_digests[header] = fragment_lengths
        
        print("\n--- 📊 Analysis Complete ---")
        max_fragments = 0
        segment_with_most = ""
        
        for header, fragments in all_digests.items():
            num_fragments = len(fragments)
            print(f"Segment '{header}': {num_fragments} fragments.")
            
            if num_fragments > max_fragments:
                max_fragments = num_fragments
                segment_with_most = header
            elif num_fragments == max_fragments:
                segment_with_most += f", {header}"

        
        print("\n--- 🏆 Result ---")
        if max_fragments == 1:
             print("No segment had more than 1 fragment (0 EcoRI sites found in any segment).")
        else:
            print(f"The segment(s) with the most fragments ({max_fragments}) is/are: '{segment_with_most}'")
        
        
        plot_combined_gel(all_digests, LADDER_SIZES)
        
        plot_separate_gels(all_digests, LADDER_SIZES)
        
        print("\nAll simulations complete.")