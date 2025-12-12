import csv
import math
import requests
from pathlib import Path
import matplotlib.pyplot as plt

# -------------------------
# CONFIG
# -------------------------
NCBI_EMAIL = "bontasadrian03@gmail.com"
WINDOW = 30

# same calibration as Ex1/Ex2
S_TEST = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
TARGET_IC = 27.53

# 10 influenza accessions (commonly used influenza A segments)
INFLUENZA_ACCESSIONS = [
    "NC_002017.1",
    "NC_002018.1",
    "NC_002019.1",
    "NC_002020.1",
    "NC_002021.1",
    "NC_002022.1",
    "NC_002023.1",
    "NC_002024.1",
    "NC_002025.1",
    "NC_002026.1",
]


# -------------------------
# NCBI download
# -------------------------
def fetch_ncbi_fasta(accession: str) -> str:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "email": NCBI_EMAIL,
    }
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    lines = r.text.splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return seq.upper().replace("N", "")


# -------------------------
# Sliding windows
# -------------------------
def sliding_windows(seq: str, w: int):
    seq = seq.upper().replace(" ", "").replace("\n", "")
    if w <= 0 or w > len(seq):
        return []
    return [seq[i:i+w] for i in range(len(seq) - w + 1)]


# -------------------------
# CG%
# -------------------------
def cg_percent(seq: str) -> float:
    if not seq:
        return 0.0
    cg = sum(1 for b in seq if b in "CG")
    return round(100.0 * cg / len(seq), 2)


# -------------------------
# Kappa IC raw + calibrated
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
    wins = sliding_windows(seq, w)
    xs = [cg_percent(win) for win in wins]
    ys = [kappa_ic(win) for win in wins]
    return xs, ys


def center_of_weight(xs, ys):
    if not xs:
        return 0.0, 0.0
    return round(sum(xs)/len(xs), 2), round(sum(ys)/len(ys), 2)


# -------------------------
# Main
# -------------------------
def main():
    here = Path(__file__).resolve().parent  # Project_L10/L10
    out_dir = here  # save next to ex1.jpg / ex1.1.jpg

    all_points_x = []
    all_points_y = []
    centers = []

    for acc in INFLUENZA_ACCESSIONS:
        print(f"[*] Downloading {acc} ...")
        seq = fetch_ncbi_fasta(acc)
        print(f"[+] {acc} length = {len(seq)}")

        if len(seq) < WINDOW:
            print(f"[-] Skip {acc} (too short)")
            continue

        xs, ys = pattern(seq, WINDOW)
        cx, cy = center_of_weight(xs, ys)

        all_points_x.extend(xs)
        all_points_y.extend(ys)
        centers.append((acc, cx, cy, len(seq), len(xs)))

    # Save centers CSV
    centers_csv = out_dir / "ex3_influenza_centers.csv"
    with centers_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["accession", "center_CG", "center_IC", "seq_len", "num_windows"])
        for row in centers:
            w.writerow(row)
    print(f"[+] Saved {centers_csv.name}")

    # Plot 1: stain
    plt.figure(figsize=(10, 6))
    plt.scatter(all_points_x, all_points_y, s=8)
    plt.title("Digital stains of influenza genomes (all windows)")
    plt.xlabel("(C+G)%")
    plt.ylabel("Kappa IC (calibrated)")
    plt.grid(True)
    plt.tight_layout()
    stain_png = out_dir / "ex3_influenza_stain.png"
    plt.savefig(stain_png, dpi=200)
    print(f"[+] Saved {stain_png.name}")

    # Plot 2: centers
    plt.figure(figsize=(10, 6))
    xs_cent = [c[1] for c in centers]
    ys_cent = [c[2] for c in centers]
    labels = [c[0] for c in centers]

    plt.scatter(xs_cent, ys_cent)
    for x, y, lab in zip(xs_cent, ys_cent, labels):
        plt.annotate(lab, (x, y))

    plt.title("Centers of digital stains (influenza)")
    plt.xlabel("Center (C+G)%")
    plt.ylabel("Center (Kappa IC)")
    plt.grid(True)
    plt.tight_layout()
    centers_png = out_dir / "ex3_influenza_centers.png"
    plt.savefig(centers_png, dpi=200)
    print(f"[+] Saved {centers_png.name}")

    plt.show()


if __name__ == "__main__":
    main()
