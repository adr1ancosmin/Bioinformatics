# ex2.py
import sys
from collections import Counter
import matplotlib.pyplot as plt

# --- Genetic code (RNA codons -> amino acid 3-letter code) ---
genetic_code = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

# ---------- Helpers ----------
def parse_fasta(path):
    """Yield (header, sequence) for each entry in a FASTA file."""
    header, seq = None, []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            else:
                if header is None:
                    continue
                seq.append(line)
    if header:
        yield header, "".join(seq)

def read_concat_sequence(path):
    """Concatenate all sequences from a FASTA into a single string."""
    pieces = []
    for _, s in parse_fasta(path):
        pieces.append(s)
    return "".join(pieces)

def count_codons(sequence):
    """
    Count non-overlapping codons in frame 0.
    Skip triplets containing letters outside A/U/C/G (T -> U first).
    """
    rna = sequence.upper().replace("T", "U")
    valid = {"A", "U", "C", "G"}
    counts = Counter()
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        if set(codon) <= valid:
            counts[codon] += 1
    return counts

def top_n(counter, n=10):
    return counter.most_common(n)

def label_with_aa(codon):
    aa = genetic_code.get(codon, "?")
    return f"{codon} ({aa})"

def plot_top_n_codons(counts, n, title, outfile):
    top = top_n(counts, n)
    if not top:
        print(f"[WARN] No codons to plot for {title}")
        return
    labels = [label_with_aa(c) for c, _ in top]
    values = [v for _, v in top]

    plt.figure()
    plt.bar(labels, values)
    for i, v in enumerate(values):
        plt.text(i, v, str(v), ha="center", va="bottom")
    plt.title(title)
    plt.xlabel("Codon (Amino acid)")
    plt.ylabel("Count")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    try:
        plt.show()
    except Exception:
        pass

def plot_compare_union_topn(counts_a, counts_b, n, name_a, name_b, outfile):
    """
    Compare counts for the union of each genome's Top-n codons.
    Creates a grouped bar chart with up to 2n codons on X.
    """
    aN = [c for c, _ in top_n(counts_a, n)]
    bN = [c for c, _ in top_n(counts_b, n)]
    codons = []
    for c in aN + bN:
        if c not in codons:
            codons.append(c)
    if not codons:
        print("[WARN] Nothing to compare.")
        return

    xa = list(range(len(codons)))
    width = 0.45

    vals_a = [counts_a.get(c, 0) for c in codons]
    vals_b = [counts_b.get(c, 0) for c in codons]

    plt.figure()
    plt.bar([x - width/2 for x in xa], vals_a, width=width, label=name_a)
    plt.bar([x + width/2 for x in xa], vals_b, width=width, label=name_b)
    plt.xticks(xa, [label_with_aa(c) for c in codons], rotation=30, ha="right")
    for x, v in zip([x - width/2 for x in xa], vals_a):
        plt.text(x, v, str(v), ha="center", va="bottom")
    for x, v in zip([x + width/2 for x in xa], vals_b):
        plt.text(x, v, str(v), ha="center", va="bottom")
    plt.legend()
    plt.title(f"Comparison: Union of Top-{n} Codons")
    plt.xlabel("Codon (Amino acid)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    try:
        plt.show()
    except Exception:
        pass

def map_codons_to_amino_acids(codons):
    """Return the amino acids (with possible duplicates) for given codons."""
    return [genetic_code.get(c, "?") for c in codons]

def unique_excluding_stop(aas):
    """Order-preserving unique AAs, excluding 'Stop'."""
    out, seen = [], set()
    for aa in aas:
        if aa == "Stop":
            continue
        if aa not in seen:
            seen.add(aa)
            out.append(aa)
    return out

# ---------- Main ----------
def main():
    covid_path     = "covid19.fasta"
    influenza_path = "influenzab.fasta"
    N = 10  # Top-N codons

    if len(sys.argv) >= 3:
        covid_path, influenza_path = sys.argv[1], sys.argv[2]
    elif len(sys.argv) == 2:
        covid_path = sys.argv[1]

    # Load & count
    covid_seq = read_concat_sequence(covid_path)
    flu_seq   = read_concat_sequence(influenza_path)

    covid_counts = count_codons(covid_seq)
    flu_counts   = count_codons(flu_seq)

    # Charts: Top-N codons (one each)
    plot_top_n_codons(covid_counts, n=N,
                      title=f"COVID-19: Top {N} most frequent codons",
                      outfile=f"top{N}_covid19.png")
    plot_top_n_codons(flu_counts, n=N,
                      title=f"Influenza: Top {N} most frequent codons",
                      outfile=f"top{N}_influenza.png")

    # One comparison chart (union of both Top-N)
    plot_compare_union_topn(covid_counts, flu_counts, N,
                            name_a="COVID-19", name_b="Influenza",
                            outfile=f"compare_top{N}.png")

    # Console: Top-N codons with AA labels (Stop may appear here at codon level)
    covid_topN_codons = [c for c, _ in top_n(covid_counts, N)]
    flu_topN_codons   = [c for c, _ in top_n(flu_counts, N)]
    covid_topN_aas    = map_codons_to_amino_acids(covid_topN_codons)
    flu_topN_aas      = map_codons_to_amino_acids(flu_topN_codons)

    print(f"\n=== Top-{N} codons & amino acids ===")
    print("COVID-19:")
    for (codon, cnt), aa in zip(top_n(covid_counts, N), covid_topN_aas):
        print(f"  {codon} ({aa}) : {cnt}")
    print("Influenza:")
    for (codon, cnt), aa in zip(top_n(flu_counts, N), flu_topN_aas):
        print(f"  {codon} ({aa}) : {cnt}")

    # AI prompts: built from Top-N amino acids but EXCLUDING 'Stop'
    covid_prompt_aas = unique_excluding_stop(covid_topN_aas)
    flu_prompt_aas   = unique_excluding_stop(flu_topN_aas)

    covid_prompt = (
        "Given common foods, identify which foods are naturally very low in or "
        f"lack the amino acids: {', '.join(covid_prompt_aas)}. "
        "Return a concise table (Food | Amino acid(s) lacking | Note). "
        "Assume typical, unfortified foods and standard portion sizes."
    )
    flu_prompt = (
        "Given common foods, identify which foods are naturally very low in or "
        f"lack the amino acids: {', '.join(flu_prompt_aas)}. "
        "Return a concise table (Food | Amino acid(s) lacking | Note). "
        "Assume typical, unfortified foods and standard portion sizes."
    )

    print("\n=== AI prompts (STOP excluded) ===")
    print("COVID-19 amino-acid prompt:")
    print(covid_prompt)
    print("\nInfluenza amino-acid prompt:")
    print(flu_prompt)

if __name__ == "__main__":
    main()
