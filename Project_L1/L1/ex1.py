sequence = "ATTTCGCCGATA"

def find_alphabet(string):
    seen = {}
    for char in string:
        if char not in seen:
            seen[char] = True
    return list(seen.keys())

alphabet = find_alphabet(sequence)
print(alphabet)
