
import tkinter as tk
import os

def read_fasta(file_path):
    sequence = ""
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                sequence += line.strip().upper().replace(' ', '')
        return sequence
    except FileNotFoundError:
        return None

def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n = len(seq1)
    m = len(seq2)
    matrix = [[0] * (m + 1) for _ in range(n + 1)]
    max_score = 0
    max_pos = None

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            
            diag = matrix[i - 1][j - 1] + score
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            
            # Local alignment (Smith-Waterman) includes 0 as a choice
            matrix[i][j] = max(0, diag, up, left)
            
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)

    # Simplified traceback: We only care if a high-scoring local alignment exists
    # If max_score is above a threshold, a significant hit is found
    threshold = 10 
    
    if max_score > threshold:
        # Returns the strength of the hit and the end position
        return max_score, max_pos
    return None

def align_segments_and_draw():
    # Attempt to read files
    seq1_full = read_fasta("covid.fasta")
    seq2_full = read_fasta("influenza.fasta")

    if seq1_full is None or seq2_full is None:
        canvas.delete("all")
        canvas.create_text(250, 100, text="Error: Could not find covid.fasta or influenza.fasta in the current directory.", fill="red")
        return

    # Segmentation parameters (The "in-between layer solution")
    segment_size = 500
    step_size = 250 # 50% overlap

    n_full = len(seq1_full)
    m_full = len(seq2_full)
    
    canvas_size = 500
    canvas.delete("all")
    canvas.create_text(canvas_size/2, 10, text=f"COVID-19 vs Influenza Similarity Map (Segments: {segment_size}bp)", fill="black")
    
    # Axis labels
    canvas.create_text(canvas_size/2, canvas_size - 10, text=f"COVID-19 Genome ({n_full}bp)")
    canvas.create_text(10, canvas_size/2, text=f"Influenza Genome ({m_full}bp)", angle=90)
    
    # Alignment loop (Block-wise Smith-Waterman)
    for i_start in range(0, n_full, step_size):
        for j_start in range(0, m_full, step_size):
            i_end = min(i_start + segment_size, n_full)
            j_end = min(j_start + segment_size, m_full)
            
            segment1 = seq1_full[i_start:i_end]
            segment2 = seq2_full[j_start:j_end]
            
            # Check for a significant local alignment hit
            hit = smith_waterman(segment1, segment2)
            
            if hit:
                max_score, max_pos = hit
                
                # Normalize coordinates for the canvas (excluding margins)
                margin = 40
                draw_area = canvas_size - 2 * margin
                
                # Start and end positions of the matched block in full genome coordinates
                x1_full = i_start
                x2_full = i_end
                y1_full = j_start
                y2_full = j_end

                # Map full coordinates to canvas coordinates (X axis for Seq1, Y axis for Seq2)
                x1_canvas = margin + (x1_full / n_full) * draw_area
                x2_canvas = margin + (x2_full / n_full) * draw_area
                y1_canvas = margin + (y1_full / m_full) * draw_area
                y2_canvas = margin + (y2_full / m_full) * draw_area
                
                # Use max_score to determine color/strength
                score_norm = min(max_score / 100, 1.0) # Normalize score up to a max of 100
                red = int(255 * score_norm)
                color = f'#{red:02x}0000' # Shades of red
                
                canvas.create_rectangle(x1_canvas, y1_canvas, x2_canvas, y2_canvas, fill=color, outline='')

# Setup the main window
root = tk.Tk()
root.title("Genome Similarity Visualizer (Local Alignment)")

# Create Canvas for the visualization
canvas_size = 500
canvas = tk.Canvas(root, width=canvas_size, height=canvas_size, bg='white')
canvas.pack(padx=10, pady=10)

# Alignment button
align_button = tk.Button(root, text="Align and Visualize Genomes", command=align_segments_and_draw)
align_button.pack(pady=10)

# Initial message
canvas.create_text(canvas_size/2, canvas_size/2, text="Click 'Align and Visualize Genomes' to start.\nEnsure 'covid.fasta' and 'influenza.fasta' are in the same directory.", fill="blue")


root.mainloop()