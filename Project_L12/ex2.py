import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests

# ===============================
# CONFIG
# ===============================
NCBI_EMAIL = "bontasadrian03@gmail.com"   # NCBI recommends a real email
TOOL_NAME = "influenza_motif_scan"
DB = "nuccore"
N_GENOMES = 10
QUERY = 'Influenza A virus[Organism] AND "complete genome"[Title]'
RESULTS_DIR = "results"
REQUEST_DELAY_SEC = 0.35  # be nice to NCBI

alphabet = "ACGT"
bg = 0.25
pseudo = 1

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

# ===============================
# PWM + LOG-LIK MODEL (Task 1)
# ===============================
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

weights = (counts + pseudo) / (N + 4 * pseudo)
loglik = np.log(weights / bg)  # ln(weight/0.25)
df_loglik = pd.DataFrame(loglik, index=list(alphabet), columns=[f"pos{j+1}" for j in range(L)])

# ===============================
# NCBI HELPERS (E-utilities)
# ===============================
def ncbi_esearch(query, retmax=10):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": DB,
        "term": query,
        "retmax": retmax,
        "retmode": "json",
        "tool": TOOL_NAME,
        "email": NCBI_EMAIL,
    }
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    return data["esearchresult"]["idlist"]

def ncbi_efetch_fasta(id_):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": DB,
        "id": id_,
        "rettype": "fasta",
        "retmode": "text",
        "tool": TOOL_NAME,
        "email": NCBI_EMAIL,
    }
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    return r.text

def parse_fasta(fasta_text):
    lines = [ln.strip() for ln in fasta_text.splitlines() if ln.strip()]
    header = lines[0]
    seq = "".join(lines[1:]).upper()
    # accession = first token after ">"
    acc = header[1:].split()[0]
    return acc, header, seq

# ===============================
# SCANNING
# ===============================
def score_window(window, df_loglik):
    # if ambiguous letters exist, return NaN so we skip them
    for ch in window:
        if ch not in "ACGT":
            return np.nan
    s = 0.0
    for j, ch in enumerate(window):
        s += df_loglik.loc[ch, f"pos{j+1}"]
    return s

def scan_sequence(seq, df_loglik, L):
    scores = np.full(len(seq) - L + 1, np.nan, dtype=float)
    for i in range(len(seq) - L + 1):
        w = seq[i:i+L]
        scores[i] = score_window(w, df_loglik)
    return scores

def plot_signal(acc, scores, L, out_png):
    x = np.arange(1, len(scores) + 1)  # 1-based window start

    # best hit ignoring NaNs
    valid_idx = np.where(~np.isnan(scores))[0]
    best_i = valid_idx[np.argmax(scores[valid_idx])]
    best_start = best_i + 1
    best_score = scores[best_i]

    plt.figure()
    plt.plot(x, scores)
    plt.title(f"Influenza genome scan (acc={acc}) | motif length={L}")
    plt.xlabel("Window start position (1-based)")
    plt.ylabel("Log-likelihood score")

    # mark best location
    plt.axvline(best_start, linestyle="--")
    plt.text(best_start, best_score, f" best @ {best_start}", rotation=90, va="bottom")

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    return best_start, best_score

# ===============================
# MAIN
# ===============================
def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print(f"Searching NCBI: {QUERY}")
    ids = ncbi_esearch(QUERY, retmax=N_GENOMES)
    print(f"Found {len(ids)} IDs. Downloading {len(ids)} records...")

    summary_rows = []

    for k, id_ in enumerate(ids, start=1):
        time.sleep(REQUEST_DELAY_SEC)

        fasta = ncbi_efetch_fasta(id_)
        acc, header, seq = parse_fasta(fasta)

        print(f"[{k}/{len(ids)}] {acc} | length={len(seq)}")

        scores = scan_sequence(seq, df_loglik, L)

        # save scores CSV
        df_out = pd.DataFrame({
            "start_pos_1based": np.arange(1, len(scores) + 1),
            "score": scores
        })
        csv_path = os.path.join(RESULTS_DIR, f"{acc}_scores.csv")
        df_out.to_csv(csv_path, index=False)

        # plot
        png_path = os.path.join(RESULTS_DIR, f"{acc}_signal.png")
        best_start, best_score = plot_signal(acc, scores, L, png_path)

        summary_rows.append({
            "accession": acc,
            "ncbi_id": id_,
            "length": len(seq),
            "best_start_1based": int(best_start),
            "best_score": float(best_score),
            "scores_csv": csv_path,
            "plot_png": png_path,
            "header": header
        })

    summary = pd.DataFrame(summary_rows).sort_values("best_score", ascending=False)
    summary_path = os.path.join(RESULTS_DIR, "summary_top_hits.csv")
    summary.to_csv(summary_path, index=False)

    print("\nDONE ✅")
    print(f"Saved plots + scores in: {RESULTS_DIR}/")
    print(f"Summary file: {summary_path}")
    print("\nTop hits:")
    print(summary[["accession", "best_start_1based", "best_score", "length"]].head(10).to_string(index=False))

if __name__ == "__main__":
    main()
