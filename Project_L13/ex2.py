import json
import random
from collections import defaultdict


def generate_random_dna(length=50):
    """
    Generate a random DNA sequence.
    
    Parameters:
    -----------
    length : int
        Length of the DNA sequence (default: 50)
    
    Returns:
    --------
    str : Random DNA sequence
    """
    nucleotides = ['A', 'T', 'G', 'C']
    return ''.join(random.choice(nucleotides) for _ in range(length))


def compute_transition_probabilities(dna_sequence):
    """
    Compute transition probabilities between nucleotides in a DNA sequence.
    This creates a "mini chatGPT" style model that learns which letter follows which.
    
    Parameters:
    -----------
    dna_sequence : str
        DNA sequence to analyze
    
    Returns:
    --------
    dict : Transition probability matrix as nested dictionary
    """
    nucleotides = ['A', 'T', 'G', 'C']
    
    # Initialize transition counts
    transition_counts = defaultdict(lambda: defaultdict(int))
    
    # Count transitions (what letter follows what letter)
    for i in range(len(dna_sequence) - 1):
        current = dna_sequence[i]
        next_char = dna_sequence[i + 1]
        transition_counts[current][next_char] += 1
    
    # Convert counts to probabilities
    transition_matrix = {}
    for nucleotide in nucleotides:
        total = sum(transition_counts[nucleotide].values())
        
        if total > 0:
            # Calculate probabilities
            transition_matrix[nucleotide] = {
                target: transition_counts[nucleotide][target] / total
                for target in nucleotides
            }
        else:
            # If nucleotide never appears, use uniform distribution
            transition_matrix[nucleotide] = {
                target: 0.25 for target in nucleotides
            }
    
    return transition_matrix, transition_counts


def display_transition_matrix(transition_matrix, transition_counts):
    """
    Display the transition matrix in a readable format.
    
    Parameters:
    -----------
    transition_matrix : dict
        Transition probability matrix
    transition_counts : dict
        Raw transition counts
    """
    nucleotides = ['A', 'T', 'G', 'C']
    
    print("\n" + "="*60)
    print("TRANSITION COUNTS (Raw frequencies)")
    print("="*60)
    print("From -> To:  A    T    G    C")
    print("-" * 60)
    
    for from_nuc in nucleotides:
        counts = [transition_counts[from_nuc][to_nuc] for to_nuc in nucleotides]
        total = sum(counts)
        print(f"   {from_nuc}    ->  {counts[0]:3d}  {counts[1]:3d}  {counts[2]:3d}  {counts[3]:3d}  (total: {total})")
    
    print("\n" + "="*60)
    print("TRANSITION PROBABILITIES (Mini-ChatGPT Model)")
    print("="*60)
    print("From -> To:    A      T      G      C")
    print("-" * 60)
    
    for from_nuc in nucleotides:
        probs = [transition_matrix[from_nuc][to_nuc] for to_nuc in nucleotides]
        print(f"   {from_nuc}    ->  {probs[0]:.3f}  {probs[1]:.3f}  {probs[2]:.3f}  {probs[3]:.3f}")
    
    print("="*60)


def save_transition_matrix_json(transition_matrix, filename="transition_matrix.json"):
    """
    Save the transition matrix to a JSON file.
    
    Parameters:
    -----------
    transition_matrix : dict
        Transition probability matrix
    filename : str
        Output filename
    """
    with open(filename, 'w') as f:
        json.dump(transition_matrix, f, indent=2)
    
    print(f"\n✓ Transition matrix saved to: {filename}")


def predict_next_nucleotide(current_nucleotide, transition_matrix):
    """
    Predict the next nucleotide based on transition probabilities.
    This is how a "mini chatGPT" would generate the next character.
    
    Parameters:
    -----------
    current_nucleotide : str
        Current nucleotide
    transition_matrix : dict
        Transition probability matrix
    
    Returns:
    --------
    str : Predicted next nucleotide
    """
    nucleotides = ['A', 'T', 'G', 'C']
    probabilities = [transition_matrix[current_nucleotide][nuc] for nuc in nucleotides]
    return random.choices(nucleotides, weights=probabilities, k=1)[0]


def generate_sequence(transition_matrix, start_nucleotide, length=20):
    """
    Generate a new DNA sequence using the learned transition probabilities.
    This demonstrates the "mini chatGPT" in action.
    
    Parameters:
    -----------
    transition_matrix : dict
        Transition probability matrix
    start_nucleotide : str
        Starting nucleotide
    length : int
        Length of sequence to generate
    
    Returns:
    --------
    str : Generated DNA sequence
    """
    sequence = start_nucleotide
    current = start_nucleotide
    
    for _ in range(length - 1):
        next_nuc = predict_next_nucleotide(current, transition_matrix)
        sequence += next_nuc
        current = next_nuc
    
    return sequence


def main():
    print("="*60)
    print("DNA TRANSITION PROBABILITY CALCULATOR")
    print("Building a 'Mini-ChatGPT' for DNA sequences")
    print("="*60)
    
    # Generate random DNA sequence
    dna_sequence = generate_random_dna(length=50)
    
    print(f"\nRandom DNA Sequence ({len(dna_sequence)} nucleotides):")
    print(f"{dna_sequence}")
    
    # Compute transition probabilities
    transition_matrix, transition_counts = compute_transition_probabilities(dna_sequence)
    
    # Display results
    display_transition_matrix(transition_matrix, transition_counts)
    
    # Save to JSON
    save_transition_matrix_json(transition_matrix, "transition_matrix.json")
    
    # Demonstrate the model by generating new sequences
    print("\n" + "="*60)
    print("GENERATING NEW SEQUENCES (Mini-ChatGPT in action)")
    print("="*60)
    
    for start in ['A', 'T', 'G', 'C']:
        generated = generate_sequence(transition_matrix, start, length=20)
        print(f"Starting with '{start}': {generated}")
    
    print("\n" + "="*60)


if __name__ == "__main__":
    main()
