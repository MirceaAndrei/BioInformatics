
import tkinter as tk
from tkinter import scrolledtext

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=0):
    n = len(seq1)
    m = len(seq2)
    matrix = [[0] * (m + 1) for _ in range(n + 1)]

    # 1. Initialize the matrix
    for i in range(n + 1):
        matrix[i][0] = i * gap
    for j in range(m + 1):
        matrix[0][j] = j * gap

    # 2. Fill the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            
            diag = matrix[i - 1][j - 1] + score
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            
            matrix[i][j] = max(diag, up, left)

    # 3. Traceback (from bottom right)
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = n, m
    matches = 0
    
    while i > 0 or j > 0:
        current_score = matrix[i][j]
        
        # Check Diagonal move (Match/Mismatch)
        if i > 0 and j > 0:
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            if current_score == matrix[i - 1][j - 1] + score:
                aligned_seq1 = seq1[i - 1] + aligned_seq1
                aligned_seq2 = seq2[j - 1] + aligned_seq2
                if seq1[i - 1] == seq2[j - 1]:
                    matches += 1
                i -= 1
                j -= 1
                continue

        # Check Up move (Gap in seq2)
        if i > 0 and current_score == matrix[i - 1][j] + gap:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        
        # Check Left move (Gap in seq1)
        elif j > 0 and current_score == matrix[i][j - 1] + gap:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    alignment_length = len(aligned_seq1)
    similarity = round((matches / alignment_length) * 100) if alignment_length > 0 else 0
    
    output = f"Alignment:\n{aligned_seq1}\n{aligned_seq2}\n\n"
    output += f"Matches = {matches}\n"
    # CORRECTED LINE: Used 'alignment_length' instead of 'length'
    output += f"Length = {alignment_length}\n" 
    output += f"Similarity = {similarity} %"
    
    return output

def align_sequences():
    seq1 = seq1_entry.get().upper()
    seq2 = seq2_entry.get().upper()
    
    try:
        match = int(match_entry.get())
        mismatch = int(mismatch_entry.get())
        gap = int(gap_entry.get())
    except ValueError:
        result_text.delete('1.0', tk.END)
        result_text.insert(tk.END, "Error: Scoring parameters must be integers.")
        return

    if not seq1 or not seq2:
        result_text.delete('1.0', tk.END)
        result_text.insert(tk.END, "Error: Sequences cannot be empty.")
        return

    result = needleman_wunsch(seq1, seq2, match, mismatch, gap)
    
    result_text.delete('1.0', tk.END)
    result_text.insert(tk.END, result)

# Setup the main window
root = tk.Tk()
root.title("Needleman-Wunsch DNA Aligner")

# Create frames for organization
sequence_frame = tk.Frame(root, padx=10, pady=10)
sequence_frame.pack(fill='x')

param_frame = tk.Frame(root, padx=10, pady=10)
param_frame.pack(fill='x')

button_frame = tk.Frame(root, padx=10, pady=10)
button_frame.pack(fill='x')

result_frame = tk.Frame(root, padx=10, pady=10)
result_frame.pack(fill='both', expand=True)

# Sequence 1 Input
tk.Label(sequence_frame, text="Sequence 1 (Sq 1):").grid(row=0, column=0, sticky='w')
seq1_entry = tk.Entry(sequence_frame, width=30)
seq1_entry.insert(0, "ACCGTGAAGCCAATAC")
seq1_entry.grid(row=0, column=1, padx=5, pady=5)

# Sequence 2 Input
tk.Label(sequence_frame, text="Sequence 2 (Sq 2):").grid(row=1, column=0, sticky='w')
seq2_entry = tk.Entry(sequence_frame, width=30)
seq2_entry.insert(0, "AGCGTGCAGCCAATAC")
seq2_entry.grid(row=1, column=1, padx=5, pady=5)

# Parameters Input
tk.Label(param_frame, text="Gap:").grid(row=0, column=0, sticky='w')
gap_entry = tk.Entry(param_frame, width=5)
gap_entry.insert(0, "0")
gap_entry.grid(row=0, column=1, padx=5, pady=5)

tk.Label(param_frame, text="Mach:").grid(row=1, column=0, sticky='w')
match_entry = tk.Entry(param_frame, width=5)
match_entry.insert(0, "1")
match_entry.grid(row=1, column=1, padx=5, pady=5)

tk.Label(param_frame, text="MMach:").grid(row=2, column=0, sticky='w')
mismatch_entry = tk.Entry(param_frame, width=5)
mismatch_entry.insert(0, "-1")
mismatch_entry.grid(row=2, column=1, padx=5, pady=5)

# Align Button
align_button = tk.Button(button_frame, text="Align", command=align_sequences, width=15)
align_button.pack()

# Result Output (scrolled text area)
tk.Label(result_frame, text="Show Alignment:").pack(anchor='w')
result_text = scrolledtext.ScrolledText(result_frame, wrap=tk.WORD, width=40, height=10)
result_text.pack(fill='both', expand=True)

root.mainloop()