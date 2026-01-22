import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------
# 0) INPUT: known motif sequences
# -------------------------------
motifs = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT",
]

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

alphabet = "ACGT"
bg = 0.25  # null model: uniform A,C,G,T (like your slide)
pseudo = 1 # Laplace pseudocount (for Weight + Log-likelihood)

# -------------------------------
# 1) Count matrix
# -------------------------------
def count_matrix(motifs, alphabet="ACGT"):
    L = len(motifs[0])
    counts = np.zeros((len(alphabet), L), dtype=int)
    for m in motifs:
        if len(m) != L:
            raise ValueError("All motifs must have the same length.")
        for j, ch in enumerate(m):
            counts[alphabet.index(ch), j] += 1
    return counts

counts = count_matrix(motifs, alphabet)
L = counts.shape[1]
N = len(motifs)

df_counts = pd.DataFrame(counts, index=list(alphabet), columns=[f"pos{j+1}" for j in range(L)])

# -------------------------------
# 2) Weight matrix (with pseudocounts)
#    weight_ij = (count_ij + pseudo) / (N + 4*pseudo)
# -------------------------------
weights = (counts + pseudo) / (N + 4*pseudo)
df_weights = pd.DataFrame(weights, index=list(alphabet), columns=df_counts.columns)

# -------------------------------
# 3) Relative frequencies matrix (raw)
#    freq_ij = count_ij / N
# -------------------------------
freqs = counts / N
df_freqs = pd.DataFrame(freqs, index=list(alphabet), columns=df_counts.columns)

# -------------------------------
# 4) Log-likelihood matrix
#    ll_ij = ln( weight_ij / bg )
# -------------------------------
loglik = np.log(weights / bg)
df_loglik = pd.DataFrame(loglik, index=list(alphabet), columns=df_counts.columns)

# Pretty printing
pd.set_option("display.precision", 4)

print("\n1) COUNT MATRIX")
print(df_counts)

print("\n2) WEIGHT MATRIX (Laplace pseudocount = 1)")
print(df_weights)

print("\n3) RELATIVE FREQUENCIES MATRIX (counts/N)")
print(df_freqs)

print("\n4) LOG-LIKELIHOODS MATRIX ln(weight/0.25)")
print(df_loglik)

# -------------------------------
# 5) Scan sequence S with sliding windows of motif length
# -------------------------------
def score_window(window, df_loglik, alphabet="ACGT"):
    score = 0.0
    for j, ch in enumerate(window):
        score += df_loglik.loc[ch, f"pos{j+1}"]
    return score

scores = []
for i in range(len(S) - L + 1):
    w = S[i:i+L]
    sc = score_window(w, df_loglik, alphabet)
    scores.append((i+1, w, sc))

df_scores = pd.DataFrame(scores, columns=["start_pos_1based", "window", "score"])
df_scores_sorted = df_scores.sort_values("score", ascending=False)

print("\nTOP 10 WINDOWS BY SCORE")
print(df_scores_sorted.head(10).to_string(index=False))

# Basic "signal" check: z-score vs all windows
mean = df_scores["score"].mean()
std = df_scores["score"].std(ddof=0)
df_scores["z"] = (df_scores["score"] - mean) / std

best = df_scores_sorted.iloc[0]
print("\nSIGNAL CHECK")
print(f"Mean score = {mean:.4f}, Std = {std:.4f}")
print(f"Best window starts at {int(best['start_pos_1based'])}, window={best['window']}, score={best['score']:.4f}")

# Plot scores across positions
plt.figure()
plt.plot(df_scores["start_pos_1based"], df_scores["score"])
plt.title("Sliding-window log-likelihood score across S")
plt.xlabel("Window start position (1-based)")
plt.ylabel("Score")
plt.show()
