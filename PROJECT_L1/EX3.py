import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
from collections import Counter
import string
import time # Added for simulating a longer task

# --- Core Biological Sequence Algorithms (Integrating lab1_1 & lab1_2 logic) ---

def read_fasta(filepath: str) -> tuple[str, str]:
    """
    Reads a FASTA file, extracts the header, and concatenates the sequence lines 
    into a single string (the 'buffer').
    
    Args:
        filepath: The path to the FASTA file.

    Returns:
        A tuple (header_line, sequence_buffer). Both are empty strings on error or if no sequence is found.
    """
    header = ""
    sequence_lines = []
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    if not header:
                        header = line[1:]
                    elif sequence_lines:
                        break
                elif header:
                    sequence_lines.append(line)
        
        sequence_buffer = "".join(sequence_lines)
        return header, sequence_buffer
        
    except Exception as e:
        messagebox.showerror("File Error", f"An error occurred while reading the file: {e}")
        return "", ""


def analyze_sequence(sequence: str, progress_bar: ttk.Progressbar) -> str:
    """
    Applies the letter analysis algorithms to the sequence buffer,
    updating the progress bar during the process.
    
    Args:
        sequence: The concatenated biological sequence string (the buffer).
        progress_bar: The ttk.Progressbar widget to update.
        
    Returns:
        A formatted string with the analysis results.
    """
    progress_bar.configure(value=0, mode='determinate')
    progress_bar.start()

    if not sequence:
        progress_bar.stop()
        return "No biological sequence found to analyze."

    # Pre-process the sequence
    cleaned_sequence = "".join(c.upper() for c in sequence if c.isalpha())
    total_letters = len(cleaned_sequence)
    
    if total_letters == 0:
        progress_bar.stop()
        return "Sequence found, but it contained no valid alphabetical characters for analysis."
    
    # Step 1: Find unique alphabet
    unique_alphabet = sorted(list(set(cleaned_sequence)))
    alphabet_str = "".join(unique_alphabet)
    
    # Simulate partial progress for a longer process
    time.sleep(1) 
    progress_bar.configure(value=50)
    progress_bar.update()

    # Step 2: Calculate letter frequency
    counts = Counter(cleaned_sequence)

    results = [
        "--- FASTA Sequence Analysis ---",
        f"Total Letters Analyzed: {total_letters}",
        "\n[1] Sequence Alphabet (Unique Letters):",
        f"{alphabet_str}",
        "\n[2] Letter Frequency (as percentage of total letters):"
    ]
    
    for letter in sorted(counts.keys()):
        count = counts[letter]
        percentage = (count * 100.0) / total_letters
        results.append(f"  {letter} : {percentage:.2f}%")
    
    # Finalize progress
    time.sleep(1)
    progress_bar.configure(value=100)
    progress_bar.update()
    progress_bar.stop()

    return "\n".join(results)


# --- GUI Application (Using tkinter) ---

class FastaAnalyzerGUI:
    def __init__(self, master):
        self.master = master
        master.title("FASTA File Analyzer")
        master.geometry("600x500")

        # 1. Select File Button
        self.select_button = tk.Button(master, text="Choose FASTA File (*.fna / *.fasta)", 
                                       command=self.select_file, 
                                       bg="#4CAF50", fg="white", 
                                       font=("Arial", 12), padx=10, pady=5)
        self.select_button.pack(pady=20, padx=10, fill='x')

        # 2. Progress Bar
        self.progress_bar = ttk.Progressbar(master, orient='horizontal', mode='determinate', length=500)
        self.progress_bar.pack(pady=(0, 10))
        self.progress_bar.pack_forget() # Initially hide the progress bar

        # 3. Output Text Area
        self.output_label = tk.Label(master, text="Analysis Results:", font=("Arial", 10, "bold"))
        self.output_label.pack(pady=(0, 5), padx=10, anchor='w')

        self.results_text = tk.Text(master, wrap='word', height=20, font=("Courier New", 10), 
                                    bg="#f0f0f0", relief="groove")
        self.results_text.pack(pady=(0, 10), padx=10, fill='both', expand=True)
        
        self.results_text.insert(tk.END, "Press the button above to select a FASTA file and run the analysis.")

    def select_file(self):
        """Opens a file dialog, reads the file, runs the analysis, and displays the results."""
        
        file_path = filedialog.askopenfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fna *.fasta"), ("All files", "*.*")],
            title="Select FASTA Sequence File"
        )

        if file_path:
            self.results_text.delete('1.0', tk.END)
            self.results_text.insert(tk.END, f"Selected file: {os.path.basename(file_path)}\n\n")
            
            # Show the progress bar
            self.progress_bar.pack(pady=(0, 10))
            self.master.update_idletasks() # Update GUI to show the progress bar

            header, sequence_buffer = read_fasta(file_path)

            if sequence_buffer:
                self.results_text.insert(tk.END, f"Header/Info Line: {header}\n")
                self.results_text.insert(tk.END, f"Raw Sequence Length (in buffer): {len(sequence_buffer)}\n\n")
                
                # Run the integrated analysis, passing the progress bar
                analysis_output = analyze_sequence(sequence_buffer, self.progress_bar)
                
                self.results_text.insert(tk.END, analysis_output)
            else:
                self.results_text.insert(tk.END, "Analysis failed. Could not find a valid sequence in the file.")
                if header:
                     self.results_text.insert(tk.END, f"\n(Header found: {header})")
            
            # Hide the progress bar after analysis is complete
            self.progress_bar.pack_forget()

if __name__ == "__main__":
    root = tk.Tk()
    app = FastaAnalyzerGUI(root)
    root.mainloop()