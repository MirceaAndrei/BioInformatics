
import tkinter as tk
from tkinter import scrolledtext

def calculate_alignment_metrics(aligned_seq1, aligned_seq2, matches, mismatch_count, gap_count, match_score, mismatch_score, gap_penalty, seq1_original_len, seq2_original_len):
    alignment_length = len(aligned_seq1)
    
    # 1. Alignment Score (Total Score)
    # Equation: Score = (Matches * Match_Score) + (Mismatches * Mismatch_Score) + (Gaps * Gap_Penalty)
    # Note: Gap_Penalty is applied per gap character in the alignment.
    # The Gap_Count already accounts for the total number of gap characters.
    alignment_score = (matches * match_score) + (mismatch_count * mismatch_score) + (gap_count * gap_penalty)

    # 2. Percent Identity
    # Equation: % Identity = (Identities / Alignment Length) * 100
    percent_identity = round((matches / alignment_length) * 100) if alignment_length > 0 else 0
    
    # 3. Percent Similarity (using a simple variant)
    # Equation: % Similarity = (Identities / Length of Shorter Original Sequence) * 100
    # This metric ignores gaps and penalizes longer alignments by comparing matches to the smallest sequence length.
    shorter_original_len = min(seq1_original_len, seq2_original_len)
    percent_similarity = round((matches / shorter_original_len) * 100) if shorter_original_len > 0 else 0
    
    return alignment_score, percent_identity, percent_similarity

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=0):
    n = len(seq1)
    m = len(seq2)
    matrix = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        matrix[i][0] = i * gap
    for j in range(m + 1):
        matrix[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            
            diag = matrix[i - 1][j - 1] + score
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            
            matrix[i][j] = max(diag, up, left)

    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = n, m
    matches = 0
    mismatch_count = 0
    gap_count = 0
    
    while i > 0 or j > 0:
        current_score = matrix[i][j]
        
        if i > 0 and j > 0:
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            if current_score == matrix[i - 1][j - 1] + score:
                aligned_seq1 = seq1[i - 1] + aligned_seq1
                aligned_seq2 = seq2[j - 1] + aligned_seq2
                if seq1[i - 1] == seq2[j - 1]:
                    matches += 1
                else:
                    mismatch_count += 1
                i -= 1
                j -= 1
                continue

        if i > 0 and current_score == matrix[i - 1][j] + gap:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            gap_count += 1
            i -= 1
        elif j > 0 and current_score == matrix[i][j - 1] + gap:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            gap_count += 1
            j -= 1

    alignment_length = len(aligned_seq1)
    
    alignment_score, percent_identity, percent_similarity = calculate_alignment_metrics(
        aligned_seq1, aligned_seq2, matches, mismatch_count, gap_count, match, mismatch, gap, n, m
    )

    output = f"Alignment:\n{aligned_seq1}\n{aligned_seq2}\n\n"
    output += f"--- Scoring Results ---\n"
    output += f"1. Alignment Score (Total Score): {alignment_score}\n"
    output += f"2. Percent Identity: {percent_identity} %\n"
    output += f"3. Percent Similarity (vs Shorter Seq): {percent_similarity} %\n"
    output += f"-----------------------\n"
    output += f"Matches = {matches}\n"
    output += f"Length = {alignment_length}\n"
    
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

    # Call the alignment function
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