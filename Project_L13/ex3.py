import json
import random
import re
from collections import defaultdict


def generate_random_text():
    """
    Generate a random English text of approximately 300 characters.
    
    Returns:
    --------
    str : Random English text
    """
    # Sample sentences to create a coherent text
    sentences = [
        "The cat sat on the mat.",
        "A dog ran through the park.",
        "The sun shines bright in the sky.",
        "Birds fly high above the trees.",
        "The moon glows softly at night.",
        "Fish swim in the deep blue ocean.",
        "Children play games in the garden.",
        "The wind blows gently through the leaves.",
        "Rain falls from the dark clouds.",
        "Stars twinkle in the evening sky.",
        "A river flows to the sea.",
        "The flowers bloom in spring.",
        "Snow covers the mountain peaks.",
        "Bees buzz around the flowers.",
        "The forest is full of life.",
    ]
    
    # Build text to approximately 300 characters
    text = ""
    while len(text) < 300:
        text += random.choice(sentences) + " "
    
    return text.strip()[:300]


def tokenize_text(text):
    """
    Tokenize text into words (lowercase, keeping punctuation separate).
    
    Parameters:
    -----------
    text : str
        Input text
    
    Returns:
    --------
    list : List of tokens (words)
    """
    # Convert to lowercase and split into words, keeping punctuation
    tokens = re.findall(r'\b\w+\b|[.,!?;]', text.lower())
    return tokens


def create_word_to_symbol_mapping(words):
    """
    Create a mapping from words to ASCII symbols for easier representation.
    
    Parameters:
    -----------
    words : list
        List of unique words
    
    Returns:
    --------
    dict, dict : word_to_symbol and symbol_to_word mappings
    """
    unique_words = sorted(set(words))
    
    # Use printable ASCII characters starting from '!'
    # ASCII printable range: 33-126
    word_to_symbol = {}
    symbol_to_word = {}
    
    for i, word in enumerate(unique_words):
        if i < 94:  # We have 94 printable ASCII characters
            symbol = chr(33 + i)  # Start from '!'
        else:
            # If we have more than 94 words, use combinations
            symbol = f"#{i}"
        
        word_to_symbol[word] = symbol
        symbol_to_word[symbol] = word
    
    return word_to_symbol, symbol_to_word


def compute_word_transition_probabilities(tokens):
    """
    Compute transition probabilities between words (Mini-GPT model).
    This learns which word typically follows which word.
    
    Parameters:
    -----------
    tokens : list
        List of tokens (words)
    
    Returns:
    --------
    dict : Transition probability matrix as nested dictionary
    dict : Raw transition counts
    """
    # Count transitions (which word follows which word)
    transition_counts = defaultdict(lambda: defaultdict(int))
    
    for i in range(len(tokens) - 1):
        current_word = tokens[i]
        next_word = tokens[i + 1]
        transition_counts[current_word][next_word] += 1
    
    # Convert counts to probabilities
    transition_matrix = {}
    
    for word in transition_counts:
        total = sum(transition_counts[word].values())
        
        if total > 0:
            transition_matrix[word] = {
                next_word: count / total
                for next_word, count in transition_counts[word].items()
            }
    
    return transition_matrix, transition_counts


def display_transition_info(transition_matrix, transition_counts, word_to_symbol, symbol_to_word):
    """
    Display transition information in a readable format.
    
    Parameters:
    -----------
    transition_matrix : dict
        Transition probability matrix
    transition_counts : dict
        Raw transition counts
    word_to_symbol : dict
        Mapping from words to symbols
    symbol_to_word : dict
        Mapping from symbols to words
    """
    print("\n" + "="*70)
    print("WORD TO SYMBOL MAPPING")
    print("="*70)
    
    # Display mapping in columns
    words = sorted(symbol_to_word.items(), key=lambda x: x[0])
    for i in range(0, len(words), 3):
        row = words[i:i+3]
        print("  ".join([f"{symbol:3s} = {word:15s}" for symbol, word in row]))
    
    print("\n" + "="*70)
    print("WORD TRANSITION COUNTS (Mini-GPT Training Data)")
    print("="*70)
    
    for word in sorted(transition_counts.keys())[:10]:  # Show first 10 for brevity
        symbol = word_to_symbol[word]
        print(f"\nWord '{word}' [{symbol}] is followed by:")
        
        # Sort by count, descending
        sorted_transitions = sorted(
            transition_counts[word].items(),
            key=lambda x: x[1],
            reverse=True
        )
        
        for next_word, count in sorted_transitions[:5]:  # Top 5
            next_symbol = word_to_symbol[next_word]
            prob = transition_matrix[word][next_word]
            print(f"  → '{next_word}' [{next_symbol}]: {count} times (prob: {prob:.3f})")
    
    if len(transition_counts) > 10:
        print(f"\n... and {len(transition_counts) - 10} more words")
    
    print("\n" + "="*70)


def save_transition_matrix_json(transition_matrix, word_to_symbol, symbol_to_word, 
                                 filename="word_transition_matrix.json"):
    """
    Save the transition matrix and mappings to a JSON file.
    
    Parameters:
    -----------
    transition_matrix : dict
        Transition probability matrix
    word_to_symbol : dict
        Word to symbol mapping
    symbol_to_word : dict
        Symbol to word mapping
    filename : str
        Output filename
    """
    output = {
        "word_to_symbol": word_to_symbol,
        "symbol_to_word": symbol_to_word,
        "transition_matrix": transition_matrix,
        "transition_matrix_with_symbols": {}
    }
    
    # Also create a version using symbols
    for word, transitions in transition_matrix.items():
        word_symbol = word_to_symbol[word]
        output["transition_matrix_with_symbols"][word_symbol] = {
            word_to_symbol[next_word]: prob
            for next_word, prob in transitions.items()
        }
    
    with open(filename, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n✓ Word transition matrix saved to: {filename}")


def generate_text(transition_matrix, start_word, max_words=20):
    """
    Generate new text using the learned word transition probabilities.
    This is the Mini-GPT in action!
    
    Parameters:
    -----------
    transition_matrix : dict
        Transition probability matrix
    start_word : str
        Starting word
    max_words : int
        Maximum number of words to generate
    
    Returns:
    --------
    str : Generated text
    """
    if start_word not in transition_matrix:
        return f"{start_word} (cannot continue - word not in training data)"
    
    words = [start_word]
    current_word = start_word
    
    for _ in range(max_words - 1):
        if current_word not in transition_matrix:
            break
        
        # Get possible next words and their probabilities
        next_words = list(transition_matrix[current_word].keys())
        probabilities = list(transition_matrix[current_word].values())
        
        # Choose next word based on probabilities
        next_word = random.choices(next_words, weights=probabilities, k=1)[0]
        words.append(next_word)
        current_word = next_word
    
    # Reconstruct text with proper spacing
    text = ""
    for word in words:
        if word in ['.', ',', '!', '?', ';']:
            text = text.rstrip() + word + " "
        else:
            text += word + " "
    
    return text.strip()


def main():
    print("="*70)
    print("WORD TRANSITION PROBABILITY CALCULATOR")
    print("Building a 'Mini-GPT' for English text")
    print("="*70)
    
    # Generate random text
    text = generate_random_text()
    
    print(f"\nRandom English Text ({len(text)} characters):")
    print("-" * 70)
    print(text)
    print("-" * 70)
    
    # Tokenize text into words
    tokens = tokenize_text(text)
    print(f"\nTokenized into {len(tokens)} words")
    print(f"Unique words: {len(set(tokens))}")
    
    # Create word-to-symbol mapping
    word_to_symbol, symbol_to_word = create_word_to_symbol_mapping(tokens)
    
    # Compute transition probabilities
    transition_matrix, transition_counts = compute_word_transition_probabilities(tokens)
    
    # Display results
    display_transition_info(transition_matrix, transition_counts, word_to_symbol, symbol_to_word)
    
    # Save to JSON
    save_transition_matrix_json(transition_matrix, word_to_symbol, symbol_to_word,
                                "word_transition_matrix.json")
    
    # Demonstrate text generation (Mini-GPT in action)
    print("\n" + "="*70)
    print("GENERATING NEW TEXT (Mini-GPT in action)")
    print("="*70)
    
    # Find common starting words
    start_words = ['the', 'a', 'in', 'and']
    available_starts = [w for w in start_words if w in transition_matrix]
    
    if not available_starts:
        available_starts = list(transition_matrix.keys())[:3]
    
    for start_word in available_starts[:3]:
        generated = generate_text(transition_matrix, start_word, max_words=15)
        print(f"\nStarting with '{start_word}':")
        print(f"  {generated}")
    
    print("\n" + "="*70)


if __name__ == "__main__":
    main()
