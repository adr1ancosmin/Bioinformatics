S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"

MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_PENALTY = 0


def needleman_wunsch(s1, s2, match=MATCH_SCORE,
                     mismatch=MISMATCH_SCORE, gap=GAP_PENALTY):
    n = len(s1)
    m = len(s2)

    # score matrix (n+1) x (m+1)
    score = [[0] * (m + 1) for _ in range(n + 1)]

    # initialize first row / column (global alignment)
    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap

    # fill the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            score[i][j] = max(diag, up, left)

    # traceback
    i, j = n, m
    align1 = []
    align2 = []

    while i > 0 or j > 0:
        if i > 0 and j > 0 and \
           score[i][j] == score[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch):
            align1.append(s1[i - 1])
            align2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and score[i][j] == score[i - 1][j] + gap:
            align1.append(s1[i - 1])
            align2.append("-")
            i -= 1
        else:
            align1.append("-")
            align2.append(s2[j - 1])
            j -= 1

    align1 = "".join(reversed(align1))
    align2 = "".join(reversed(align2))

    return align1, align2, score


def print_alignment(a1, a2):
    # middle line with '|' for matches, ' ' for mismatches/gaps
    mid = []
    for x, y in zip(a1, a2):
        mid.append("|" if x == y else " ")
    mid = "".join(mid)

    matches = sum(1 for x, y in zip(a1, a2) if x == y)
    length = len(a1)
    similarity = matches / length * 100

    print(a1)
    print(mid)
    print(a2)
    print()
    print(f"Matches   = {matches}")
    print(f"Length    = {length}")
    print(f"Similarity = {similarity:.0f} %")


if __name__ == "__main__":
    alignment1, alignment2, M = needleman_wunsch(S1, S2)
    print_alignment(alignment1, alignment2)
