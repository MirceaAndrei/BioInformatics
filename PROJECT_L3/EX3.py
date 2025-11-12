from tkinter import filedialog, Tk, Label, Button, Entry, messagebox
import matplotlib.pyplot as plt
import math

def basic_tm(S):
    A = S.count('A')
    T = S.count('T')
    G = S.count('G')
    C = S.count('C')
    return 4 * (G + C) + 2 * (A + T)

def advanced_tm(S):
    na = 0.001
    length = len(S)
    if length == 0:
        return None
    gc_count = S.count('G') + S.count('C')
    gc_perc = (gc_count / length) * 100
    log_na = math.log10(na)
    tm = 81.5 + 16.6 * log_na + 0.41 * gc_perc - (600 / length)
    return round(tm, 2)

def calculate_tm_signals(sequence, window_size=9):
    tm_basic = []
    tm_advanced = []

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        tm_basic.append(basic_tm(window))
        tm_advanced.append(advanced_tm(window))

    return tm_basic, tm_advanced

def show_main_chart(tm_basic, tm_advanced, threshold):
    x = list(range(len(tm_basic)))
    plt.figure(figsize=(10, 5))
    plt.plot(x, tm_basic, label="Wallace Formula", color='blue')
    plt.plot(x, tm_advanced, label="Nearest Neighbor", color='red')
    plt.axhline(y=threshold, color='gray', linestyle='--', label=f"Threshold = {threshold}°C")
    plt.title("Melting Temperature Profile (Window = 9)")
    plt.xlabel("Window Start Position")
    plt.ylabel("Temperature (°C)")
    plt.ylim(0, 100)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def show_threshold_chart(tm_basic, tm_advanced, threshold):
    x = list(range(len(tm_basic)))

    plt.figure(figsize=(10, 5))
    plt.fill_between(x, threshold, tm_basic, where=[val >= threshold for val in tm_basic],
                     interpolate=True, color='blue', alpha=0.5, label="Wallace > threshold")

    plt.fill_between(x, threshold, tm_advanced, where=[val >= threshold for val in tm_advanced],
                     interpolate=True, color='red', alpha=0.5, label="NN > threshold")

    plt.axhline(y=threshold, color='gray', linestyle='--', label=f"Threshold = {threshold}°C")
    plt.title("Tm Segments Above Threshold")
    plt.xlabel("Window Start Position")
    plt.ylabel("Temperature (°C)")
    plt.ylim(0, 100)
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

            threshold_input = entry_threshold.get().strip()
            try:
                threshold = float(threshold_input)
            except ValueError:
                messagebox.showerror("Error", "Invalid threshold. Please enter a number.")
                return

            tm_basic, tm_advanced =calculate_tm_signals(sequence, window_size=9)

            print(f"Wallace Formula: Min = {min(tm_basic)} °C, Max = {max(tm_basic)} °C")
            print(f"Nearest Neighbor: Min = {min(tm_advanced)} °C, Max = {max(tm_advanced)} °C")

            show_main_chart(tm_basic, tm_advanced, threshold)
            show_threshold_chart(tm_basic, tm_advanced, threshold)

    except Exception as e:
        messagebox.showerror("Error", f"Failed to process file:\n{e}")

root = Tk()
root.title("DNA Melting Temperature Analyzer")
root.geometry("500x250")

label = Label(root, text="Select a FASTA file to analyze", font=("Arial", 12))
label.pack(pady=10)

button = Button(root, text="Open FASTA File", font=("Arial", 12), command=open_file)
button.pack()

label_thresh = Label(root, text="Set Threshold (°C) for Filtered View:", font=("Arial", 11))
label_thresh.pack(pady=10)

entry_threshold = Entry(root, font=("Arial", 11))
entry_threshold.insert(0, "5")
entry_threshold.pack()

root.mainloop()
