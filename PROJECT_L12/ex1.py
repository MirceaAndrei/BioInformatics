import math

matrix = {
    'A': [0.3, 0.6, 0.1, 0.0, 0.0, 0.6, 0.7, 0.2, 0.1],
    'C': [0.2, 0.2, 0.1, 0.0, 0.0, 0.2, 0.1, 0.1, 0.2],
    'G': [0.1, 0.1, 0.7, 1.0, 0.0, 0.1, 0.1, 0.5, 0.1],
    'T': [0.4, 0.1, 0.1, 0.0, 1.0, 0.1, 0.1, 0.2, 0.6]
}

pseudo = 0.001
background = 0.25
pwm = {}

for nt, probs in matrix.items():
    pwm[nt] = []
    for p in probs:
        adjusted_p = (p + pseudo) / (1 + 4 * pseudo)
        pwm[nt].append(math.log(adjusted_p / background))

def calculate_score(sequence_window):
    return sum(pwm[nt][i] for i, nt in enumerate(sequence_window))

s = "CCGTAGGTAGGTCCG"
window_size = 9

for i in range(len(s) - window_size + 1):
    window = s[i:i+window_size]
    score = calculate_score(window)
    print(f"Pos {i}: {window} | Score: {score:.4f}")