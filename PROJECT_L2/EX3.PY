import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline

# Function to read FASTA file and extract sequence
def read_fasta(file_path):
    sequence = []
    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):  # Skip header lines
                sequence.append(line.strip().upper())
    return "".join(sequence)

# Compute relative frequencies (%) for each sliding window
def compute_frequencies(sequence, window_size=1):
    if len(sequence) < window_size:
        messagebox.showwarning("Warning", "Sequence is shorter than window size.")
        return {}

    alphabet = sorted(set(sequence))
    results = {base: [] for base in alphabet}

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        total = len(window)
        for base in alphabet:
            freq_percent = (window.count(base) / total) * 100
            results[base].append(freq_percent)

    return results

# Smooth with cubic spline interpolation
def smooth_curve(x, y, num_points=500):
    if len(x) < 4:
        return x, y  # too short to interpolate
    spline = make_interp_spline(x, y, k=3)  # cubic spline
    x_smooth = np.linspace(x.min(), x.max(), num_points)
    y_smooth = spline(x_smooth)
    return x_smooth, y_smooth

# Plot results (continuous smooth curves)
def plot_frequencies(results):
    if not results:
        return

    plt.figure(figsize=(10, 6))
    x = np.arange(len(next(iter(results.values()))))

    for base, freqs in results.items():
        x_smooth, y_smooth = smooth_curve(x, np.array(freqs))
        plt.plot(x_smooth, y_smooth, label=base, linewidth=2)

    plt.title("Sliding Window Base Percentages (Window = 30)")
    plt.xlabel("Window Position")
    plt.ylabel("Relative Frequency (%)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

# Handle file selection and analysis
def load_and_analyze():
    file_path = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")]
    )

    if not file_path:
        return

    try:
        sequence = read_fasta(file_path)
        if not sequence:
            messagebox.showerror("Error", "No sequence data found in file.")
            return

        results = compute_frequencies(sequence, window_size=30)
        plot_frequencies(results)

    except Exception as e:
        messagebox.showerror("Error", f"Failed to process file:\n{e}")

# --- GUI setup ---
root = tk.Tk()
root.title("FASTA Sliding Window Analyzer")
root.geometry("400x200")

title_label = tk.Label(root, text="FASTA Sliding Window Analyzer", font=("Arial", 14, "bold"))
title_label.pack(pady=20)

info_label = tk.Label(root, text="Select a FASTA file to analyze (Window = 30)", font=("Arial", 10))
info_label.pack(pady=10)

analyze_button = tk.Button(root, text="Choose FASTA File", command=load_and_analyze, font=("Arial", 12))
analyze_button.pack(pady=20)

root.mainloop()