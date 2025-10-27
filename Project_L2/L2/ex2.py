S = "ABAACTGTTGCTGTTGTGC"

def find_kmers(seq, k):
    kmers = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1
    return kmers

def calculate_frequency(seq, k):
    kmers = find_kmers(seq, k)
    total = sum(kmers.values())
    result = {}
    for kmer, count in kmers.items():
        result[kmer] = f"{100 * count / total:.4f}%"
    return dict(sorted(result.items()))

print(f"Sequence: {S}")
print(f"Length: {len(S)}")
print(f"Di-nucleotides: {calculate_frequency(S, 2)}")
print(f"Tri-nucleotides: {calculate_frequency(S, 3)}")
