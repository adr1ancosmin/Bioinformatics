S = "ACTGTGTCCCGTAAGCCAATCTGTTG"

nucleotides = ['A', 'T', 'G', 'C']

dinucleotides = []
for n1 in nucleotides:
    for n2 in nucleotides:
        dinucleotides.append(n1 + n2)

trinucleotides = []
for n1 in nucleotides:
    for n2 in nucleotides:
        for n3 in nucleotides:
            trinucleotides.append(n1 + n2 + n3)

def calculate_frequency(sequence, kmers):
    frequencies = {}
    total_count = 0
    
    for kmer in kmers:
        count = 0
        kmer_len = len(kmer)
        
        for i in range(len(sequence) - kmer_len + 1):
            if sequence[i:i+kmer_len] == kmer:
                count += 1
        
        frequencies[kmer] = count
        total_count += count
    
    relative_frequencies = {}
    if total_count > 0:
        for kmer, count in frequencies.items():
            relative_frequencies[kmer] = count / total_count
    else:
        relative_frequencies = {kmer: 0.0 for kmer in kmers}
    
    return frequencies, relative_frequencies

print("DI-NUCLEOTIDES")
di_counts, di_relative = calculate_frequency(S, dinucleotides)
print(f"Sequence: {S}")
print(f"Length: {len(S)}")
print(f"Total possible di-nucleotide positions: {len(S) - 1}")

for kmer in sorted(dinucleotides):
    if di_counts[kmer] > 0:
        print(f"{kmer}: count={di_counts[kmer]}, percentage={100*di_relative[kmer]:.4f}%")

print("\nTRI-NUCLEOTIDES")
tri_counts, tri_relative = calculate_frequency(S, trinucleotides)
print(f"Total possible tri-nucleotide positions: {len(S) - 2}")

for kmer in sorted(trinucleotides):
    if tri_counts[kmer] > 0:
        print(f"{kmer}: count={tri_counts[kmer]}, percentage={100*tri_relative[kmer]:.4f}%")

print("\nSUMMARY")
print(f"Total di-nucleotides found: {sum(di_counts.values())}")
print(f"Unique di-nucleotides found: {sum(1 for c in di_counts.values() if c > 0)}/{len(dinucleotides)}")
print(f"Total tri-nucleotides found: {sum(tri_counts.values())}")
print(f"Unique tri-nucleotides found: {sum(1 for c in tri_counts.values() if c > 0)}/{len(trinucleotides)}")
