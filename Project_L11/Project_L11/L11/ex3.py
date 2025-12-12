import requests
import textwrap

# ============================================================
# EXERCISE 3
# We reuse the local alignment from exercise 2 (Smith–Waterman)
# and ADD 3 different similarity scoring equations:
#
# Let:
#   M = number of matches
#   X = number of mismatches (both positions are bases)
#   G = number of positions with at least one gap ("-")
#   L_align = length of the alignment (including gaps)
#
# 1) Score1: percent identity over full alignment
#    S1 = 100 * M / L_align
#
# 2) Score2: percent identity ignoring gaps
#    S2 = 100 * M / (M + X)
#
# 3) Score3: gap-penalized similarity (gap = 0.5 mismatch)
#    S3 = 100 * M / (M + X + 0.5 * G)
# ============================================================

# ---------- CONFIG (CHANGE THESE) ----------

NCBI_EMAIL = "bontasadrian03@gmail.com"  # your email

INFLUENZA_ACC = "NC_002026.1"   # Influenza A HA segment
COVID_ACC     = "NC_045512.2"   # SARS-CoV-2 reference genome

# Smith–Waterman scoring
MATCH_SCORE    = 2
MISMATCH_SCORE = -1
GAP_PENALTY    = -2

# how big is each region (window) for step-by-step alignment
CHUNK_SIZE = 2000


# ---------- STEP 1: DOWNLOAD SEQUENCES FROM NCBI ----------

def fetch_ncbi_fasta(accession: str) -> str:
    """
    Download a nucleotide sequence from NCBI in FASTA format
    and return only the sequence (no header, no newlines).
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "email": NCBI_EMAIL,
    }

    print(f"[*] Downloading {accession} from NCBI...")
    r = requests.get(url, params=params)
    r.raise_for_status()

    lines = r.text.splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    seq = seq.upper().replace("N", "")  # remove Ns, capitalize
    print(f"[+] Downloaded {accession}, length = {len(seq)} bp")
    return seq


# ---------- STEP 2: SMITH–WATERMAN LOCAL ALIGNMENT ----------

def smith_waterman(s1: str, s2: str,
                   match: int = MATCH_SCORE,
                   mismatch: int = MISMATCH_SCORE,
                   gap: int = GAP_PENALTY):
    """
    Basic Smith–Waterman local alignment.
    Returns: aligned_s1, aligned_s2, max_score
    """

    n, m = len(s1), len(s2)
    # H[i][j] = best local score ending at s1[i-1], s2[j-1]
    H = [[0] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = H[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            delete = H[i - 1][j] + gap
            insert = H[i][j - 1] + gap
            H[i][j] = max(0, diag, delete, insert)

            if H[i][j] > max_score:
                max_score = H[i][j]
                max_pos = (i, j)

    # traceback from max_pos until we hit 0
    i, j = max_pos
    a1, a2 = [], []

    while i > 0 and j > 0 and H[i][j] != 0:
        score_current = H[i][j]
        score_diag = H[i - 1][j - 1]
        score_up = H[i - 1][j]
        score_left = H[i][j - 1]

        if score_current == score_diag + (match if s1[i - 1] == s2[j - 1] else mismatch):
            a1.append(s1[i - 1])
            a2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif score_current == score_up + gap:
            a1.append(s1[i - 1])
            a2.append("-")
            i -= 1
        elif score_current == score_left + gap:
            a1.append("-")
            a2.append(s2[j - 1])
            j -= 1
        else:
            break

    a1 = "".join(reversed(a1))
    a2 = "".join(reversed(a2))
    return a1, a2, max_score


# ---------- STEP 2.5: PRINT + 3 SCORING EQUATIONS ----------

def print_alignment(a1: str, a2: str, max_width: int = 80):
    """
    Pretty-print alignment in blocks of max_width
    and compute 3 different similarity scores.
    """

    # middle line: '|' only when both bases are equal and not gaps
    mid = "".join(
        "|" if x == y and x != "-" and y != "-" else " "
        for x, y in zip(a1, a2)
    )

    # split in blocks for nicer printing
    blocks1 = textwrap.wrap(a1, max_width)
    blocks2 = textwrap.wrap(a2, max_width)
    blocksM = textwrap.wrap(mid, max_width)

    for b1, bm, b2 in zip(blocks1, blocksM, blocks2):
        print(b1)
        print(bm)
        print(b2)
        print()

    # ----- STATS FOR SCORING -----
    matches = 0
    mismatches = 0
    gaps = 0

    for x, y in zip(a1, a2):
        if x == "-" or y == "-":
            gaps += 1
        elif x == y:
            matches += 1
        else:
            mismatches += 1

    alignment_len = len(a1)  # includes gaps

    # Score 1: percent identity over full alignment (gaps counted in length)
    if alignment_len > 0:
        score1 = 100.0 * matches / alignment_len
    else:
        score1 = 0.0

    # Score 2: percent identity ignoring gap positions
    comparable = matches + mismatches  # positions where both are bases
    if comparable > 0:
        score2 = 100.0 * matches / comparable
    else:
        score2 = 0.0

    # Score 3: gap-penalized similarity (gap = half mismatch)
    denom3 = matches + mismatches + 0.5 * gaps
    if denom3 > 0:
        score3 = 100.0 * matches / denom3
    else:
        score3 = 0.0

    print(f"Matches      = {matches}")
    print(f"Mismatches   = {mismatches}")
    print(f"Gaps         = {gaps}")
    print(f"AlignmentLen = {alignment_len}")
    print("-" * 60)
    print(f"Score1 (ID over full alignment)   = {score1:.2f}%")
    print(f"Score2 (ID without gaps)          = {score2:.2f}%")
    print(f"Score3 (gap-penalized similarity) = {score3:.2f}%")
    print("-" * 60)


# ---------- STEP 3: ALIGNMENT "LAYER" – BIG REGIONS ----------

def stepwise_alignment(seq1: str, seq2: str, chunk_size: int = CHUNK_SIZE):
    """
    Align the two genomes step by step on big regions.
    This is the "in-between layer":
    instead of one huge alignment, we align successive windows.
    """

    min_len = min(len(seq1), len(seq2))
    print(f"\n[!] Stepwise local alignment with chunk size = {chunk_size} bp")
    print(f"    Using {min_len // chunk_size + 1} regions (approx.)\n")

    region_idx = 1
    for start in range(0, min_len, chunk_size):
        end = min(start + chunk_size, min_len)
        sub1 = seq1[start:end]
        sub2 = seq2[start:end]

        print(f"=== REGION {region_idx}: bases {start+1}-{end} ===")
        a1, a2, score = smith_waterman(sub1, sub2)
        print(f"Local max score = {score}")
        print_alignment(a1, a2, max_width=80)
        region_idx += 1


# ---------- MAIN ----------

def main():
    flu = fetch_ncbi_fasta(INFLUENZA_ACC)
    covid = fetch_ncbi_fasta(COVID_ACC)

    print("\nInfluenza length:", len(flu))
    print("COVID-19 length :", len(covid))

    stepwise_alignment(flu, covid, chunk_size=CHUNK_SIZE)


if __name__ == "__main__":
    main()
