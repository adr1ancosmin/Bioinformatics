import random
import re
import sys
import ssl
import urllib.request
import urllib.error
import matplotlib.pyplot as plt


N_SAMPLES = 2000
LEN_MIN = 100
LEN_MAX = 150
OVERLAP_MIN = 10
N_PLOT = 120
GLOBAL_SEED = 7

DNA = (
"TTGACCGATCTTAGGACCGTTAACGGCTTAGGACTTACCGGAATTCGACCTTAGGCATACCGGTTAGACCTTAGGACCG"
"TTACGGAATTCGATCGGATACCTTAGGCCAATTCGATCGAACCGGATTACTTAGGACCGATTCGACCGGATACCTTACG"
"GAATTCGACCGGATTAGGCTTACCGGATTACTCGGAATTCGACCTTAGGACCGATACCGGATTAGGCTTACGGATTAC"
"TCGGAATTCGACTTAGGACCGATACCGGATTAGGCTTACCGGATTACTCGGAATTCGACTTAGGACCGATACCGGATT"
"AGGCTTACCGGATTACTCGGAATTCGACTTAGGACCGATACCGGATTAGGCTTACCGGATTACTCGGAATTCGACCTT"
"AGGACCGATACCGGATTAGGCTTACCGGATTACTCGGAATTCGACCGGATACCTTAGGACCGATTACTCGGAATTCGA"
"CTTAGGACCGATACCGGATTAGGCTTACCGGATTACTCGGAATTCGACCTTAGGACCGATACCGGATTAGGCTTACCG"
"GA"
)


def generate_reads(dna: str, n_reads: int, lmin: int, lmax: int, seed: int):
    rng = random.Random(seed)
    coords = []
    reads = []
    max_start = len(dna) - lmax
    for _ in range(n_reads):
        s = rng.randint(0, max_start)
        read_len = rng.randint(lmin, lmax)
        e = s + read_len
        reads.append(dna[s:e])
        coords.append((s, e))
    return reads, coords


def greedy_assembly(reads, min_overlap: int) -> str:
    if not reads:
        return ""
    contig = reads[0]
    for r in reads[1:]:
        lim = min(len(r), len(contig))
        joined = False
        for k in range(lim, min_overlap - 1, -1):
            if contig.endswith(r[:k]):
                contig += r[k:]
                joined = True
                break
        if not joined:
            pass
    return contig


def plot_reads(coords, genome_len: int, n_show: int, seed: int, title: str):
    rng = random.Random(seed)
    subset = coords if len(coords) <= n_show else rng.sample(coords, n_show)

    plt.figure(figsize=(14, 6))
    for idx, (s, e) in enumerate(subset):
        plt.hlines(idx, s, e, linewidth=2)
    plt.xlabel(f"Position on reference (0 → {genome_len})")
    plt.ylabel("Read index")
    plt.title(title)
    plt.tight_layout()
    plt.show()


def main():
    random.seed(GLOBAL_SEED)

    dna = DNA
    print("DNA length:", len(dna))

    reads, coords = generate_reads(
        dna, n_reads=N_SAMPLES, lmin=LEN_MIN, lmax=LEN_MAX, seed=GLOBAL_SEED
    )

    contig = greedy_assembly(reads, OVERLAP_MIN)

    print("\nORIGINAL DNA (first 200 nt):")
    print(dna[:200])
    print("\nRECONSTRUCTED DNA (first 200 nt):")
    print(contig[:200])

    print("\nOriginal length:", len(dna))
    print("Reconstructed length:", len(contig))
    print("Reconstructed is substring of original:", contig in dna)
    print("Original is substring of reconstructed:", dna in contig)

    plot_reads(
        coords,
        genome_len=len(dna),
        n_show=N_PLOT,
        seed=GLOBAL_SEED,
        title="Random fragments sampled from DNA",
    )

    summary = f"""
Assembly summary:

- Chaining uses ≥{OVERLAP_MIN} overlap; many fragments cannot merge.
- Random sequence still contains repeated short motifs causing ambiguity.
- Contig is shorter than the real (synthetic) DNA.

More advanced assemblers use graph structures and scoring rules to improve contiguity.
"""
    print(summary)


if __name__ == "__main__":
    main()
