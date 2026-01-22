import math
from collections import defaultdict

# =========================
# INPUT SEQUENCES (GIVEN)
# =========================
S1 = "ATCGATTCGATATCATACACGTAT"     # CpG island (+)
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"   # Non-island (-)
S_test = "CAGGTTGGAAACGTAA"

NUCLEOTIDES = ['A', 'C', 'G', 'T']


# =========================
# TRANSITION COUNT
# =========================
def count_transitions(sequence):
    counts = {n: defaultdict(int) for n in NUCLEOTIDES}

    for i in range(len(sequence) - 1):
        from_n = sequence[i]
        to_n = sequence[i + 1]
        counts[from_n][to_n] += 1

    return counts


# =========================
# TRANSITION PROBABILITIES
# =========================
def transition_probabilities(counts):
    probs = {n: {} for n in NUCLEOTIDES}

    for from_n in NUCLEOTIDES:
        total = sum(counts[from_n].values())
        for to_n in NUCLEOTIDES:
            if total == 0:
                probs[from_n][to_n] = 0
            else:
                probs[from_n][to_n] = counts[from_n][to_n] / total

    return probs


# =========================
# LOG-LIKELIHOOD MATRIX
# log2( P_plus / P_minus )
# =========================
def log_likelihood_matrix(p_plus, p_minus):
    ll = {n: {} for n in NUCLEOTIDES}

    for i in NUCLEOTIDES:
        for j in NUCLEOTIDES:
            if p_plus[i][j] == 0 or p_minus[i][j] == 0:
                ll[i][j] = 0
            else:
                ll[i][j] = math.log(p_plus[i][j] / p_minus[i][j], 2)

    return ll


# =========================
# SEQUENCE SCORING
# =========================
def score_sequence(sequence, ll_matrix):
    score = 0
    for i in range(len(sequence) - 1):
        score += ll_matrix[sequence[i]][sequence[i + 1]]
    return score


# =========================
# MAIN
# =========================
if __name__ == "__main__":
    # Count transitions
    counts_plus = count_transitions(S1)
    counts_minus = count_transitions(S2)

    # Probabilities
    probs_plus = transition_probabilities(counts_plus)
    probs_minus = transition_probabilities(counts_minus)

    # Log-likelihood matrix
    ll_matrix = log_likelihood_matrix(probs_plus, probs_minus)

    # Score test sequence
    final_score = score_sequence(S_test, ll_matrix)

    # Output
    print("Log-likelihood score:", round(final_score, 4))

    if final_score > 0:
        print("Result: CpG ISLAND (+)")
    else:
        print("Result: NON-CpG ISLAND (-)")
