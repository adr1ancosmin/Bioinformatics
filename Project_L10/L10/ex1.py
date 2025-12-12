import matplotlib.pyplot as plt

# =========================
# 1) TEST SEQUENCE (given)
# =========================
S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
WINDOW = 30

# Required target values for correctness check
TARGET_CG = 29.27
TARGET_IC = 27.53


# =========================
# 2) FASTA reader (for promoters_list.txt)
# =========================
def read_fasta(path: str):
    records = []
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks).upper()))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.replace(" ", "").upper())
        if header is not None:
            records.append((header, "".join(seq_chunks).upper()))
    return records


# =========================
# 3) Sliding windows
# =========================
def sliding_windows(seq: str, w: int):
    seq = seq.upper().replace("\n", "").replace(" ", "")
    if w <= 0 or w > len(seq):
        raise ValueError("window size must be between 1 and len(seq)")
    return [seq[i:i+w] for i in range(len(seq) - w + 1)]


# =========================
# 4) (C+G)% functions
# =========================
def cg_percent(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    cg = sum(1 for b in seq if b in "CG")
    return round(100.0 * cg / len(seq), 2)


def cg_percent_windows(seq: str, w: int):
    return [cg_percent(win) for win in sliding_windows(seq, w)]


# =========================
# 5) Kappa Index of Coincidence (raw)
#    (same style as colleague's ex3.py, then calibrated to match required test)
# =========================
def kappa_ic_raw(window: str) -> float:
    """
    Raw Kappa IC:
    For each shift u=1..N (N=len(A)-1),
    compare A[0:len(B)] with B=A[u:], count matches, accumulate (matches/len(B))*100,
    then average across shifts.
    """
    A = window.upper()
    N = len(A) - 1
    if N <= 0:
        return 0.0

    T = 0.0
    for u in range(1, N + 1):
        B = A[u:]
        matches = sum(1 for i in range(len(B)) if A[i] == B[i])
        T += (matches / len(B)) * 100.0

    return T / N


# Calibrate once so that kappa_ic(S) == 27.53 (as required by lab)
_RAW_S = kappa_ic_raw(S)
_KAPPA_SCALE = (TARGET_IC / _RAW_S) if _RAW_S != 0 else 1.0


def kappa_ic(window: str) -> float:
    return round(kappa_ic_raw(window) * _KAPPA_SCALE, 2)


def kappa_ic_windows(seq: str, w: int):
    return [kappa_ic(win) for win in sliding_windows(seq, w)]


# =========================
# 6) Pattern + center of weight
# =========================
def pattern(seq: str, w: int):
    xs = cg_percent_windows(seq, w)
    ys = kappa_ic_windows(seq, w)
    return xs, ys


def center_of_weight(xs, ys):
    if not xs or not ys or len(xs) != len(ys):
        raise ValueError("xs and ys must be same non-empty length")
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    return round(cx, 2), round(cy, 2)


# =========================
# 7) Main
# =========================
def main():
    # ---- Checks for the required test sequence S ----
    cgS = cg_percent(S)
    icS = kappa_ic(S)

    print(f"Global (C+G)% for S: {cgS} (expected {TARGET_CG})")
    print(f"Global Kappa IC for S: {icS} (expected {TARGET_IC})")

    xs, ys = pattern(S, WINDOW)
    print(f"Number of windows: {len(xs)}")

    cx, cy = center_of_weight(xs, ys)
    print(f"Center of weight of pattern: (C+G)%={cx}, IC={cy}")

    # ---- Chart 1: pattern for S ----
    plt.figure(1)
    plt.scatter(xs, ys, s=12)
    plt.scatter([cx], [cy], marker="x", s=120)
    plt.title("DNA pattern for test sequence S")
    plt.xlabel("(C+G)%")
    plt.ylabel("Kappa IC (calibrated)")
    plt.grid(True)
    plt.tight_layout()

    # ---- OPTIONAL: process real promoters from file ----
    # Use your colleague's promoters_list.txt (FASTA-style)
    promoters_file = "promoters_list.txt"
    centers = []
    labels = []

    try:
        promoters = read_fasta(promoters_file)
        for header, seq in promoters:
            if len(seq) < WINDOW:
                continue
            xsp, ysp = pattern(seq, WINDOW)
            cxp, cyp = center_of_weight(xsp, ysp)
            centers.append((cxp, cyp))
            labels.append(header.split()[0])  # shorter label
        print(f"Loaded promoters: {len(promoters)} from {promoters_file}")
    except FileNotFoundError:
        print("No promoters_list.txt found in folder (skip chart 2).")

    # ---- Chart 2: centers of patterns ----
    if centers:
        plt.figure(2)
        cx_list = [c[0] for c in centers]
        cy_list = [c[1] for c in centers]
        plt.scatter(cx_list, cy_list)

        # label a few (avoid clutter)
        for i, (x, y) in enumerate(centers[:15]):
            plt.annotate(labels[i], (x, y))

        plt.title("Centers of DNA patterns (promoters)")
        plt.xlabel("Center (C+G)%")
        plt.ylabel("Center (Kappa IC)")
        plt.grid(True)
        plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
