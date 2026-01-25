import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

S1 = "ATCGATTCGATATCATACACGTAT"
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"
S_test = "CAGGTTGGAAACGTAA"
plt.close('all')
def count_transitions(seq):
    t = defaultdict(int)
    total = defaultdict(int)
    for i in range(len(seq) - 1):
        a, b = seq[i], seq[i+1]
        t[(a, b)] += 1
        total[a] += 1
    return t, total

t_pos, total_pos = count_transitions(S1)
t_neg, total_neg = count_transitions(S2)

nucs = ['A', 'C', 'G', 'T']
mat_pos = pd.DataFrame(index=nucs, columns=nucs, dtype=float)
mat_neg = pd.DataFrame(index=nucs, columns=nucs, dtype=float)

for a in nucs:
    for b in nucs:
        mat_pos.loc[a, b] = t_pos[(a, b)] / total_pos[a] if total_pos[a] > 0 else 0
        mat_neg.loc[a, b] = t_neg[(a, b)] / total_neg[a] if total_neg[a] > 0 else 0

epsilon = 1e-10
LLM = np.log2((mat_pos + epsilon) / (mat_neg + epsilon))

def score_sequence(seq, log_matrix):
    score = 0.0
    for i in range(len(seq) - 1):
        a, b = seq[i], seq[i+1]
        score += log_matrix.loc[a, b]
    return score

score = score_sequence(S_test, LLM)

print("Observed matrix (+):")
print(mat_pos.round(2))
print("\nExpected matrix (-):")
print(mat_neg.round(2))
print("\nLog-likelihood matrix (LLM):")
print(LLM.round(2))
print("\nSequence:", S_test)
print("Score:", round(score, 4))
if score > 0:
    print("The sequence is likely to be from a CpG island.")
else:
    print("The sequence is likely NOT from a CpG island.")

plt.figure(figsize=(8,1))
plt.text(0.1, 0.5, r'$LLM = \log_2 \left( \frac{p_{i+1,j+1}}{b_{i+1,j+1}} \right)$', fontsize=18)
plt.axis('off')
plt.show()

plt.figure(figsize=(6,5))
sns.heatmap(mat_pos, annot=True, cmap="YlGnBu", fmt=".2f", cbar=True)
plt.title("Observed Matrix (+)")
plt.show()

plt.figure(figsize=(6,5))
sns.heatmap(mat_neg, annot=True, cmap="YlOrRd", fmt=".2f", cbar=True)
plt.title("Expected Matrix (-)")
plt.show()

plt.figure(figsize=(6,5))
sns.heatmap(LLM, annot=True, cmap="coolwarm", fmt=".2f", cbar=True)
plt.title("Log-Likelihood Matrix (LLM)")
plt.show()