import sys

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

def translate_gene(rna_sequence):
    amino_acid_sequence = []
    
    rna = rna_sequence.upper().strip().replace("T", "U")
    
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        
        if len(codon) < 3:
            break 
            
        amino_acid = genetic_code.get(codon) 
        
        if amino_acid:
            if amino_acid == 'Stop':
                break
            amino_acid_sequence.append(amino_acid)
        else:
            amino_acid_sequence.append('X') 
            
    return "-".join(amino_acid_sequence)

def parse_fasta(filename):
    """
    Parses a FASTA file and yields (header, sequence) tuples.
    Handles multi-line sequences.
    """
    header = None
    sequence_lines = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: 
                    continue
                
                if line.startswith('>'):
                    if header:
                        yield (header, "".join(sequence_lines))
                    
                    header = line[1:]   
                    sequence_lines = []
                else:
                    if header: 
                        sequence_lines.append(line)
            
            if header:
                yield (header, "".join(sequence_lines))
                
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.", file=sys.stderr)
        return

if __name__ == "__main__":
    fasta_filename = "sequence.fasta"
    
    print(f"Translating all sequences from '{fasta_filename}'...\n")
    print("=" * 50)
    
    sequences_found = False
    for header, sequence in parse_fasta(fasta_filename):
        sequences_found = True
        
        print(f"> {header}")
        
        try:
            protein = translate_gene(sequence)
            
            if protein:
                print(f"Protein: {protein}")
            else:
                print("Protein: (No translatable codons found)")
                
        except Exception as e:
            print(f"Could not translate sequence: {e}")
            
        print("-" * 50)

    if not sequences_found:
        print(f"No valid FASTA sequences were found in '{fasta_filename}'.")
        print("Please make sure the file exists and is in the correct format.")