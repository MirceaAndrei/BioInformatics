from tkinter import filedialog, Tk, Label, Button, messagebox
import matplotlib.pyplot as plt
import math

def basic_tm(S):
    A = S.count('A')
    T = S.count('T')
    G = S.count('G')
    C = S.count('C')
    return 4 * (G + C) + 2 * (A + T)

def advanced_tm(S, na=0.05):
    length = len(S)
    if length == 0:
        return None
    gc_count = S.count('G') + S.count('C')
    gc_perc = (gc_count / length) * 100
    log_na = math.log10(na)
    tm = 81.5 + 16.6 * log_na + 0.41 * gc_perc - (600 / length)
    return round(tm, 2)

def calculate_tm_signals(sequence, window_size=9, na=0.001):
    tm_basic = []
    tm_advanced = []

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        tm_basic.append(basic_tm(window))
        tm_advanced.append(advanced_tm(window, na))

    return tm_basic, tm_advanced

def show_chart(tm_basic, tm_advanced):
    plt.plot(tm_basic, label="Basic Tm (4GC+2AT)")
    plt.plot(tm_advanced, label="Advanced Tm (Salt-adjusted)")
    plt.title("Melting Temperature (Tm) - Sliding Window = 9")
    plt.xlabel("Window Start Position")
    plt.ylabel("Temperature (°C)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def open_file():
    filepath = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa *.txt")])
    if not filepath:
        return
    try:
        with open(filepath, "r") as file:
            lines = file.readlines()
            sequence = "".join([line.strip().upper() for line in lines if not line.startswith(">")])
            sequence = "".join([base for base in sequence if base in ['A', 'T', 'C', 'G']])

            if len(sequence) < 9:
                messagebox.showerror("Error", "Sequence too short. Minimum 9 bases required.")
                return

            tm_basic, tm_advanced = calculate_tm_signals(sequence, window_size=9, na=0.001)
            show_chart(tm_basic, tm_advanced)

    except Exception as e:
        messagebox.showerror("Error", f"Failed to process file:\n{e}")

# GUI setup
root = Tk()
root.title("Melting Temperature (Tm) Analyzer")
root.geometry("400x200")

label = Label(root, text="Select a FASTA file to analyze", font=("Arial", 12))
label.pack(pady=20)

button = Button(root, text="Open FASTA File", font=("Arial", 12), command=open_file)
button.pack()

root.mainloop()
