
import numpy as np
import matplotlib.pyplot as plt

S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
WINDOW_LENGTH = 30

def calculate_cg_percent(sequence):
    N = len(sequence)
    if N == 0:
        return 0.0
    C_count = sequence.count('C')
    G_count = sequence.count('G')
    return ((C_count + G_count) / N) * 100.0

def calculate_ic(sequence):
    N = len(sequence)
    if N < 2:
        return 0.0
    bases = ['A', 'C', 'G', 'T']
    counts = {base: sequence.count(base) for base in bases}
    
    ic_sum = 0
    for count in counts.values():
        ic_sum += count * (count - 1)
        
    ic_standard = ic_sum / (N * (N - 1))
    return ic_standard * 100.0

def run_sliding_window_analysis(sequence, window):
    cg_values = []
    ic_values = []
    window_centers = []
    
    for i in range(len(sequence) - window + 1):
        window_sequence = sequence[i:i + window]
        
        cg = calculate_cg_percent(window_sequence)
        ic = calculate_ic(window_sequence)
        
        cg_values.append(cg)
        ic_values.append(ic)
        window_centers.append(i + window / 2)
        
    return np.array(window_centers), np.array(cg_values), np.array(ic_values)

def calculate_center_of_weight(positions, values):
    numerator = np.sum(positions * values)
    denominator = np.sum(values)
    if denominator == 0:
        return 0
    return numerator / denominator

def plot_pattern(positions, cg_values, ic_values):
    fig, ax1 = plt.subplots(figsize=(12, 6))

    ax1.set_xlabel('Window Center Position (bp)')
    ax1.set_ylabel('C+G %', color='blue')
    ax1.plot(positions, cg_values, color='blue', label='C+G %')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.grid(axis='y', linestyle='--', alpha=0.7)

    ax2 = ax1.twinx()
    ax2.set_ylabel('Index of Coincidence (Scaled)', color='red')
    ax2.plot(positions, ic_values, color='red', linestyle='--', label='Index of Coincidence')
    ax2.tick_params(axis='y', labelcolor='red')

    fig.suptitle('DNA Promoter Pattern: C+G % and Index of Coincidence')
    plt.legend(loc='upper right')
    plt.show()
    return fig

def plot_center_of_weight(positions, cg_values, cow):
    plt.figure(figsize=(12, 3))
    plt.plot(positions, cg_values, color='green', alpha=0.5, label='C+G % Pattern')
    plt.axvline(x=cow, color='red', linestyle='-', linewidth=2, label=f'Center of Weight: {cow:.2f}')
    plt.scatter([cow], [np.max(cg_values) / 2], color='red', marker='X', s=200, zorder=5, label='CoW Point')
    plt.xlabel('Window Center Position (bp)')
    plt.ylabel('C+G %')
    plt.title('Center of Weight of the C+G % Pattern')
    plt.legend()
    plt.grid(True)
    plt.show()

window_centers, cg_results, ic_results = run_sliding_window_analysis(S, WINDOW_LENGTH)

avg_cg = np.mean(cg_results)
avg_ic = np.mean(ic_results)

cow_cg = calculate_center_of_weight(window_centers, cg_results)

print("--- Analysis Summary ---")
print(f"DNA Sequence Length: {len(S)} bp")
print(f"Sliding Window Length: {WINDOW_LENGTH} bp")
print(f"Number of Windows Analyzed: {len(window_centers)}")
print("\n--- Required Test Values (Overall Average of Window Results) ---")
print(f"Calculated C+G % Average: {avg_cg:.2f}")
print(f"Required C+G % Target: 29.27")
print(f"Calculated Kappa Index of Coincidence Average (Standard IC x 100): {avg_ic:.2f}")
print(f"Required IC Target: 27.53")
print("\n--- Center of Weight Calculation ---")
print(f"Center of Weight of the C+G % Pattern: {cow_cg:.2f} bp")

plot_pattern(window_centers, cg_results, ic_results)
plot_center_of_weight(window_centers, cg_results, cow_cg)