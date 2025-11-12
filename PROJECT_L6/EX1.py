import random
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Reusing the get_random_fragment_lengths function from before
def get_random_fragment_lengths(seq_length, num_fragments, min_len, max_len_user):
    """Generates a list of random fragment lengths."""
    fragment_lengths = []
    
    max_len_possible = min(seq_length, max_len_user)

    if min_len > max_len_possible:
        print(f"Error: Minimum fragment length ({min_len}) is greater than "
              f"the max possible length ({max_len_possible}).")
        return []

    print(f"Generating {num_fragments} random fragment lengths between {min_len} bp "
          f"and {max_len_possible} bp...")

    for _ in range(num_fragments):
        frag_len = random.randint(min_len, max_len_possible)
        fragment_lengths.append(frag_len)
        
    return fragment_lengths

def plot_gel_with_pyplot_single_sample_lane(sample_lengths, ladder_sizes):
    """
    Simulates the gel electrophoresis with a stylized, band-like representation,
    plotting all sample fragments in a single lane next to the ladder.
    """
    
    print("\n--- ⚡ Simulating Stylized Gel with All Samples in One Lane ⚡ ---")
    
    fig, ax = plt.subplots(figsize=(6, 10)) # Adjusted figure size

    # --- Configuration ---
    LANE_WIDTH = 0.8 # Width of each band in x-direction
    BAND_HEIGHT = 0.03 # Relative height of the band (visual thickness)
    LANE_SPACING = 1.0 # Spacing between lane centers
    
    # Now there are only 2 lanes: Ladder and Combined Samples
    NUM_LANES = 2 

    # Determine Y-axis range based on all fragment lengths for consistent scaling
    all_lengths = sample_lengths + ladder_sizes
    min_len = min(all_lengths)
    max_len = max(all_lengths)
    
    # Use a common normalization scale for all fragments
    log_min_all = math.log10(min_len)
    log_max_all = math.log10(max_len)

    def get_normalized_y_pos(length):
        if log_max_all == log_min_all: # Handle case of all same size
            return 0.5
        # Scales from 0 (shortest) to 1 (longest)
        return (math.log10(length) - log_min_all) / (log_max_all - log_min_all)
    
    # --- Set up the Gel background ---
    ax.set_facecolor('black') 
    ax.set_xlim(0, NUM_LANES * LANE_SPACING + 0.5) 
    ax.set_ylim(-0.1, 1.1) 
    
    # Hide axis ticks and labels for a cleaner look
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False) 

    # --- Draw Wells (top of the gel) ---
    for i in range(NUM_LANES):
        x_center = (i * LANE_SPACING) + LANE_SPACING/2
        well = patches.Rectangle((x_center - LANE_WIDTH/2, 1.05), 
                                 LANE_WIDTH, 0.05, 
                                 facecolor='gray', edgecolor='white', lw=1.5, clip_on=False)
        ax.add_patch(well)
    ax.text(NUM_LANES * LANE_SPACING / 2, 1.15, 'Wells (Negative Electrode)', 
            ha='center', va='bottom', fontsize=10, color='gray')


    # --- Draw Ladder Lane (first lane) ---
    x_ladder_center = LANE_SPACING / 2
    for size in ladder_sizes:
        y_pos = get_normalized_y_pos(size)
        
        ax.barh(y=y_pos, width=LANE_WIDTH, height=BAND_HEIGHT, 
                left=x_ladder_center - LANE_WIDTH/2, color='white', 
                edgecolor='darkgray', lw=0.5)
        
        ax.text(x_ladder_center - LANE_WIDTH/2 - 0.1, y_pos, 
                f'{size} bp -', 
                ha='right', va='center', color='black', fontsize=10,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, boxstyle='round,pad=0.2'))

    # --- Draw Combined Sample Lane (second lane) ---
    x_sample_center = 1 * LANE_SPACING + LANE_SPACING / 2 # The second lane
    
    # To show "more samples are located" in a region, we can draw slightly
    # thicker or brighter bands if multiple fragments fall very close.
    # For now, let's just plot all fragments in the same lane.
    
    # Collect y_positions for samples and sort them for easier processing if needed
    sample_y_positions = sorted([get_normalized_y_pos(length) for length in sample_lengths])

    # Plot each sample band
    for y_pos in sample_y_positions:
        ax.barh(y=y_pos, width=LANE_WIDTH, height=BAND_HEIGHT,
                left=x_sample_center - LANE_WIDTH/2, color='white',
                edgecolor='darkgray', lw=0.5)
                
    # To indicate regions where more samples are located,
    # we could potentially use a histogram or density plot over the y-axis,
    # but for a gel-like visual, drawing individual bands is still best.
    # If bands are very close, they will visually merge or appear denser.


    # --- Add lane labels ---
    ax.text(x_ladder_center, 0.98, 'Ladder', ha='center', va='top', color='gray', fontsize=10)
    ax.text(x_sample_center, 0.98, 'Samples', ha='center', va='top', color='gray', fontsize=10)


    # --- Add general gel labels ---
    ax.text(NUM_LANES * LANE_SPACING / 2, -0.05, 'Positive Electrode', 
            ha='center', va='top', fontsize=10, color='gray')

    ax.set_title('Simulated DNA Gel Electrophoresis', fontsize=16, color='black', pad=30)
    
    plt.tight_layout()
    plt.savefig('stylized_gel_single_sample_lane.png', facecolor='white')
    
    print("Stylized gel plot with single sample lane saved as stylized_gel_single_sample_lane.png")

# --- Main Program Execution ---
if __name__ == "__main__":
    
    MAIN_SEQ_LENGTH = 1541 
    NUM_SAMPLES = 10
    MIN_LENGTH = 100
    MAX_LENGTH = 3000
    
    LADDER_SIZES = [3000, 2000, 1500, 1200, 1000, 800, 600, 500, 400, 300, 200, 100]

    fragment_lengths = get_random_fragment_lengths(
        MAIN_SEQ_LENGTH, NUM_SAMPLES, MIN_LENGTH, MAX_LENGTH
    )
        
    if fragment_lengths:
        print("Generated fragment lengths (bp):", fragment_lengths)
        
        plot_gel_with_pyplot_single_sample_lane(fragment_lengths, LADDER_SIZES)