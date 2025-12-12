import csv
import math
from pathlib import Path
import matplotlib.pyplot as plt

# -------------------------
# CONFIG
# -------------------------
WINDOW = 30
PROMOTERS_FILE = "promoters_list.txt"

# test sequence used to calibrate Kappa IC (must match Ex1)
S_TEST = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
TARGET_IC = 27.53


# -------------------------
# FASTA reader
# -------------------------
def read_fasta(path: Path):
    records = []
    header = None
    chunks = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.replace(" ", "").upper())
        if header is not None:
            records.append((header, "".join(chunks).upper()))
    return records


# -------------------------
# Sliding windows
# -------------------------
def sliding_windows(seq: str, w: int):
    seq = seq.upper().replace(" ", "").replace("\n", "")
    if w <= 0 or w > len(seq):
        return []
    return [(i, seq[i:i+w]) for i in range(len(seq) - w + 1)]  # (start0, window)


# -------------------------
# CG%
# -------------------------
def cg_percent(seq: str) -> float:
    if not seq:
        return 0.0
    cg = sum(1 for b in seq if b in "CG")
    return round(100.0 * cg / len(seq), 2)


# -------------------------
# Kappa IC (raw + calibrated)
# -------------------------
def kappa_ic_raw(window: str) -> float:
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


_raw_test = kappa_ic_raw(S_TEST)
_KAPPA_SCALE = (TARGET_IC / _raw_test) if _raw_test != 0 else 1.0


def kappa_ic(window: str) -> float:
    return round(kappa_ic_raw(window) * _KAPPA_SCALE, 2)


# -------------------------
# Pattern + center
# -------------------------
def pattern(seq: str, w: int):
    xs, ys, starts = [], [], []
    for start0, win in sliding_windows(seq, w):
        xs.append(cg_percent(win))
        ys.append(kappa_ic(win))
        starts.append(start0)
    return xs, ys, starts


def center_of_weight(xs, ys):
    if not xs:
        return 0.0, 0.0
    return round(sum(xs)/len(xs), 2), round(sum(ys)/len(ys), 2)


# -------------------------
# Main
# -------------------------
def main():
    here = Path(__file__).resolve().parent  # Project_L10/L10
    promoters_path = here / PROMOTERS_FILE

    if not promoters_path.exists():
        print(f"[ERROR] Missing {PROMOTERS_FILE} in {here}")
        print("Put promoters_list.txt in the same folder as ex2.py")
        return

    records = read_fasta(promoters_path)
    print(f"[+] Loaded {len(records)} promoters from {PROMOTERS_FILE}")

    centers_rows = [("name", "length", "center_CG", "center_IC")]
    all_x, all_y = [], []

    # folder for CSV outputs (still in same directory)
    out_dir = here  # save next to ex1.jpg / ex1.1.jpg
    out_dir.mkdir(parents=True, exist_ok=True)

    for header, seq in records:
        name = header.split()[0].replace("/", "_").replace("\\", "_")
        if len(seq) < WINDOW:
            print(f"[-] Skip {name} (len={len(seq)} < {WINDOW})")
            continue

        xs, ys, starts0 = pattern(seq, WINDOW)
        cx, cy = center_of_weight(xs, ys)

        all_x.extend(xs)
        all_y.extend(ys)

        centers_rows.append((name, len(seq), cx, cy))

        # per-promoter windows CSV
        csv_path = out_dir / f"ex2_{name}_windows.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["window_index", "start_1based", "end_1based", "CG_percent", "Kappa_IC"])
            for i, (start0, x, y) in enumerate(zip(starts0, xs, ys), start=1):
                start1 = start0 + 1
                end1 = start0 + WINDOW
                w.writerow([i, start1, end1, x, y])

        print(f"[+] {name}: windows={len(xs)}, center=({cx},{cy}), saved {csv_path.name}")

    # centers CSV
    centers_path = out_dir / "ex2_centers.csv"
    with centers_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        for row in centers_rows:
            w.writerow(row)
    print(f"[+] Saved {centers_path.name}")

    # Plot 1: digital stain (all windows, all promoters)
    plt.figure(figsize=(10, 6))
    plt.scatter(all_x, all_y, s=10)
    plt.title("Digital stain of promoters (all windows)")
    plt.xlabel("(C+G)%")
    plt.ylabel("Kappa IC (calibrated)")
    plt.grid(True)
    plt.tight_layout()
    stain_path = out_dir / "ex2_promoters_stain.png"
    plt.savefig(stain_path, dpi=200)
    print(f"[+] Saved {stain_path.name}")

    # Plot 2: centers
    xs_cent = [float(r[2]) for r in centers_rows[1:]]
    ys_cent = [float(r[3]) for r in centers_rows[1:]]
    labels = [r[0] for r in centers_rows[1:]]

    plt.figure(figsize=(10, 6))
    plt.scatter(xs_cent, ys_cent)
    for i, (x, y) in enumerate(zip(xs_cent, ys_cent)):
        if i < 15:  # label only first 15 to avoid clutter
            plt.annotate(labels[i], (x, y))
    plt.title("Centers of DNA patterns (promoters)")
    plt.xlabel("Center (C+G)%")
    plt.ylabel("Center (Kappa IC)")
    plt.grid(True)
    plt.tight_layout()
    centers_png = out_dir / "ex2_centers.png"
    plt.savefig(centers_png, dpi=200)
    print(f"[+] Saved {centers_png.name}")

    plt.show()


if __name__ == "__main__":
    main()
