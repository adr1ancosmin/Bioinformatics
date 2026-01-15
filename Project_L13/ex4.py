import json
import random
import os


def bags_to_transition_matrix(state_bags):
    """
    Convert bag/urn representation to transition probability matrix.
    This demonstrates the teacher's concept: A = {AAAAAABBBBBA}
    
    Parameters:
    -----------
    state_bags : dict
        Dictionary where keys are states and values are bags (lists) of next states
        Example: {'A': ['A','A','A','A','A','A','B','B','B','B','B','A']}
    
    Returns:
    --------
    dict : Transition probability matrix
    """
    transition_matrix = {}
    
    for state, bag in state_bags.items():
        # Count occurrences in the bag
        total = len(bag)
        transitions = {}
        
        for next_state in set(bag):
            count = bag.count(next_state)
            transitions[next_state] = count / total
        
        transition_matrix[state] = transitions
    
    return transition_matrix


def simulate_with_bags(state_bags, initial_state, num_steps=10):
    """
    Simulate random picking from bags (teacher's ABB concept).
    At each step, randomly pick a letter from the current state's bag.
    
    Parameters:
    -----------
    state_bags : dict
        Dictionary of state bags
    initial_state : str
        Starting state
    num_steps : int
        Number of steps to simulate
    
    Returns:
    --------
    list : Sequence of states
    """
    sequence = [initial_state]
    current_state = initial_state
    
    for _ in range(num_steps):
        # Randomly pick from the bag for current state
        bag = state_bags[current_state]
        next_state = random.choice(bag)
        sequence.append(next_state)
        current_state = next_state
    
    return sequence


def display_bag_model(state_bags):
    """
    Display the bag/urn model as the teacher explained.
    Shows bags like: A = {AAAAAABBBBBA}
    """
    print("\n" + "="*70)
    print("BAG/URN MODEL (Teacher's Concept)")
    print("="*70)
    print("\nEach state has a bag with letters inside.")
    print("To get the next state, randomly pick a letter from the current bag!\n")
    
    for state, bag in state_bags.items():
        bag_str = ''.join(bag)
        print(f"{state} = {{ {bag_str} }}")
        
        # Count what's in the bag
        counts = {}
        for item in bag:
            counts[item] = counts.get(item, 0) + 1
        
        parts = [f"{count} {s}'s" if count > 1 else f"{count} {s}" for s, count in sorted(counts.items())]
        print(f"     Contains: {', '.join(parts)} (total: {len(bag)})")
        
        # Show probabilities
        for s, count in sorted(counts.items()):
            prob = count / len(bag)
            print(f"     -> {prob*100:.1f}% chance to pick '{s}'")
        print()


def demonstrate_bag_picking():
    """
    Demonstrate the teacher's bag model concept with A and B states.
    """
    print("\n" + "="*70)
    print("EXERCISE 4: SEQUENCE SYNTHESIS USING BAG MODEL")
    print("="*70)
    print("\nTeacher's example: Random picking from bags (ABB concept)")
    
    # Teacher's example bags
    state_bags = {
        'A': ['A', 'A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B', 'A'],
        'B': ['B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'A']
    }
    
    display_bag_model(state_bags)
    
    # Convert to transition matrix
    transition_matrix = bags_to_transition_matrix(state_bags)
    
    print("="*70)
    print("TRANSITION MATRIX (From Bags)")
    print("="*70)
    print("\nFrom -> To:      A        B")
    print("-" * 70)
    for state in ['A', 'B']:
        probs = [transition_matrix[state].get(s, 0.0) for s in ['A', 'B']]
        print(f"    {state}    ->   {probs[0]:.3f}    {probs[1]:.3f}")
    
    # Simulate random sequences by picking from bags
    print("\n" + "="*70)
    print("RANDOM SEQUENCES (Picking from bags)")
    print("="*70)
    print()
    
    for i in range(5):
        sequence = simulate_with_bags(state_bags, 'A', num_steps=15)
        seq_str = ''.join(sequence)
        print(f"Sequence {i+1}: {seq_str}")
        print(f"            Explanation: Start at A, randomly pick from each bag")
    
    # State vector analysis
    print("\n" + "="*70)
    print("STATE VECTOR ANALYSIS")
    print("="*70)
    
    all_states = []
    for _ in range(100):
        seq = simulate_with_bags(state_bags, 'A', num_steps=20)
        all_states.extend(seq)
    
    count_A = all_states.count('A')
    count_B = all_states.count('B')
    total = len(all_states)
    
    print(f"\nAfter many simulations (total: {total} state visits):")
    print(f"  State A: {count_A} visits ({count_A/total*100:.1f}%)")
    print(f"  State B: {count_B} visits ({count_B/total*100:.1f}%)")
    print(f"\nState vector: [A: {count_A/total:.3f}, B: {count_B/total:.3f}]")


def load_transition_matrix_json(filename):
    """
    Load a transition matrix from a JSON file.
    
    Parameters:
    -----------
    filename : str
        Name of the JSON file
    
    Returns:
    --------
    dict : Loaded transition matrix data
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File '{filename}' not found. Please run ex2.py or ex3.py first.")
    
    with open(filename, 'r') as f:
        data = json.load(f)
    
    return data


def synthesize_dna_sequence(transition_matrix, start_nucleotide='A', length=50, num_sequences=5):
    """
    Synthesize DNA sequences using the transition matrix.
    
    Parameters:
    -----------
    transition_matrix : dict
        Transition probability matrix for DNA
    start_nucleotide : str
        Starting nucleotide (default: 'A')
    length : int
        Length of sequence to generate
    num_sequences : int
        Number of sequences to generate
    
    Returns:
    --------
    list : List of generated DNA sequences
    """
    sequences = []
    
    for _ in range(num_sequences):
        sequence = start_nucleotide
        current = start_nucleotide
        
        for _ in range(length - 1):
            if current not in transition_matrix:
                # If current nucleotide has no transitions, pick random
                current = random.choice(['A', 'T', 'G', 'C'])
            else:
                # Get next nucleotide based on probabilities
                next_options = list(transition_matrix[current].keys())
                probabilities = list(transition_matrix[current].values())
                
                if sum(probabilities) > 0:
                    current = random.choices(next_options, weights=probabilities, k=1)[0]
                else:
                    current = random.choice(['A', 'T', 'G', 'C'])
            
            sequence += current
        
        sequences.append(sequence)
    
    return sequences


def synthesize_text_sequence(transition_matrix, start_word='the', max_words=30, num_sequences=5):
    """
    Synthesize text sequences using the word transition matrix.
    
    Parameters:
    -----------
    transition_matrix : dict
        Transition probability matrix for words
    start_word : str
        Starting word
    max_words : int
        Maximum number of words to generate
    num_sequences : int
        Number of sequences to generate
    
    Returns:
    --------
    list : List of generated text sequences
    """
    sequences = []
    
    for _ in range(num_sequences):
        words = [start_word]
        current_word = start_word
        
        for _ in range(max_words - 1):
            if current_word not in transition_matrix:
                break
            
            # Get possible next words and their probabilities
            next_words = list(transition_matrix[current_word].keys())
            probabilities = list(transition_matrix[current_word].values())
            
            if not next_words:
                break
            
            # Choose next word based on probabilities
            next_word = random.choices(next_words, weights=probabilities, k=1)[0]
            words.append(next_word)
            current_word = next_word
            
            # Stop at sentence end
            if next_word in ['.', '!', '?']:
                break
        
        # Format text properly
        text = ""
        for word in words:
            if word in ['.', ',', '!', '?', ';', ':']:
                text = text.rstrip() + word + " "
            else:
                text += word + " "
        
        sequences.append(text.strip())
    
    return sequences


def display_dna_sequences(sequences, start_nucleotide):
    """
    Display synthesized DNA sequences.
    
    Parameters:
    -----------
    sequences : list
        List of DNA sequences
    start_nucleotide : str
        Starting nucleotide
    """
    print("\n" + "="*70)
    print(f"SYNTHESIZED DNA SEQUENCES (starting with '{start_nucleotide}')")
    print("="*70)
    
    for i, seq in enumerate(sequences, 1):
        print(f"\nSequence {i} ({len(seq)} nucleotides):")
        # Print in chunks of 50 for readability
        for j in range(0, len(seq), 50):
            print(f"  {seq[j:j+50]}")
    
    print("\n" + "="*70)


def display_text_sequences(sequences, start_word):
    """
    Display synthesized text sequences.
    
    Parameters:
    -----------
    sequences : list
        List of text sequences
    start_word : str
        Starting word
    """
    print("\n" + "="*70)
    print(f"SYNTHESIZED TEXT SEQUENCES (starting with '{start_word}')")
    print("="*70)
    
    for i, text in enumerate(sequences, 1):
        word_count = len(text.split())
        print(f"\nSequence {i} ({word_count} words):")
        # Wrap text at 65 characters for readability
        words = text.split()
        line = "  "
        for word in words:
            if len(line) + len(word) + 1 > 67:
                print(line)
                line = "  " + word + " "
            else:
                line += word + " "
        if line.strip():
            print(line)
    
    print("\n" + "="*70)


def process_dna_matrix():
    """Process DNA transition matrix from ex2."""
    print("\n" + "="*70)
    print("DNA SEQUENCE SYNTHESIS (from Exercise 2)")
    print("="*70)
    
    try:
        data = load_transition_matrix_json("transition_matrix.json")
        transition_matrix = data
        
        print("\n✓ Loaded DNA transition matrix from 'transition_matrix.json'")
        print(f"  Matrix contains transitions for: {', '.join(transition_matrix.keys())}")
        
        # Generate sequences for each starting nucleotide
        for start_nuc in ['A', 'T', 'G', 'C']:
            sequences = synthesize_dna_sequence(
                transition_matrix, 
                start_nucleotide=start_nuc, 
                length=60, 
                num_sequences=3
            )
            display_dna_sequences(sequences, start_nuc)
        
        return True
        
    except FileNotFoundError as e:
        print(f"\n✗ Error: {e}")
        return False


def process_word_matrix():
    """Process word transition matrix from ex3."""
    print("\n" + "="*70)
    print("TEXT SEQUENCE SYNTHESIS (from Exercise 3)")
    print("="*70)
    
    try:
        data = load_transition_matrix_json("word_transition_matrix.json")
        
        # Extract the transition matrix
        if "transition_matrix" in data:
            transition_matrix = data["transition_matrix"]
            word_to_symbol = data.get("word_to_symbol", {})
        else:
            transition_matrix = data
            word_to_symbol = {}
        
        print("\n✓ Loaded word transition matrix from 'word_transition_matrix.json'")
        print(f"  Matrix contains {len(transition_matrix)} unique words")
        
        if word_to_symbol:
            print(f"  Word-to-symbol mapping available: {len(word_to_symbol)} words")
        
        # Find good starting words
        common_starts = ['the', 'a', 'an', 'in', 'on', 'at', 'to', 'for', 'of', 'and']
        available_starts = [w for w in common_starts if w in transition_matrix]
        
        if not available_starts:
            # Use first few words from the matrix
            available_starts = list(transition_matrix.keys())[:3]
        
        # Generate sequences for different starting words
        for start_word in available_starts[:4]:
            sequences = synthesize_text_sequence(
                transition_matrix,
                start_word=start_word,
                max_words=25,
                num_sequences=5
            )
            display_text_sequences(sequences, start_word)
        
        return True
        
    except FileNotFoundError as e:
        print(f"\n✗ Error: {e}")
        return False


def main():
    # Part 1: Teacher's bag model concept (ABB random picking)
    demonstrate_bag_picking()
    
    # Part 2: Synthesize from JSON transition matrices
    print("\n\n" + "="*70)
    print("PART 2: SYNTHESIS FROM JSON TRANSITION MATRICES")
    print("="*70)
    
    # Try to process DNA matrix (from ex2)
    dna_success = process_dna_matrix()
    
    # Try to process word matrix (from ex3)
    word_success = process_word_matrix()
    
    # Summary
    print("\n" + "="*70)
    print("SYNTHESIS SUMMARY")
    print("="*70)
    
    if dna_success:
        print("✓ DNA sequences synthesized successfully from JSON")
    else:
        print("✗ DNA synthesis skipped (run ex2.py first to generate transition_matrix.json)")
    
    if word_success:
        print("✓ Text sequences synthesized successfully from JSON")
    else:
        print("✗ Text synthesis skipped (run ex3.py first to generate word_transition_matrix.json)")
    
    if not dna_success and not word_success:
        print("\n⚠ Note: Run ex2.py and/or ex3.py first to generate transition matrices!")
    
    print("\n" + "="*70)
    print("COMPLETE: Bag model + JSON matrix synthesis demonstrated")
    print("="*70)


if __name__ == "__main__":
    main()
