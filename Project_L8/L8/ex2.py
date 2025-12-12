import json

IN_JSON = "te_dataset.json"


def find_all_occurrences(seq: str, pattern: str):
    """
    Return list of (start,end) 1-based inclusive for all occurrences of pattern in seq.
    """
    hits = []
    start = 0
    while True:
        idx = seq.find(pattern, start)
        if idx == -1:
            break
        s = idx + 1
        e = idx + len(pattern)
        hits.append((s, e))
        start = idx + 1  # allow overlaps if any
    return hits


def main():
    with open(IN_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)

    genome = data["genome"]
    te_library = data["te_library"]  # name -> TE full sequence

    print(f"[+] Loaded genome length: {len(genome)} bp")
    print(f"[+] TE types in library: {list(te_library.keys())}\n")

    found = []
    for name, te_seq in te_library.items():
        hits = find_all_occurrences(genome, te_seq)
        if not hits:
            print(f"[-] {name}: NOT FOUND")
        else:
            for (s, e) in hits:
                found.append((name, s, e, e - s + 1))
                print(f"[+] {name}: start={s}, end={e}, len={e - s + 1}")

    # optional: compare with ground truth if present
    if "tes_inserted" in data:
        print("\nGround truth:")
        for p in data["tes_inserted"]:
            print(f"  {p['name']}: {p['start']}-{p['end']} (len={p['length']})")


if __name__ == "__main__":
    main()
