sequence = "ATTTCGCCGATA"

def find_alphabet(string):
    seen = {}
    for char in string:
        if char not in seen:
            seen[char] = True
    return list(seen.keys())

alphabet = find_alphabet(sequence)

def relative_freq(alphabet, string):
    result = {}
    for c in alphabet:
        result[c] = str(100 * len([n for n in string if(n == c)])/len(string)) + "%"
    return result

print(relative_freq(alphabet=alphabet, string=sequence))
