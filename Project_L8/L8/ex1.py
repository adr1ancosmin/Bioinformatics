import random
import json

# ---------------- CONFIG ----------------
OUT_JSON = "te_dataset.json"
OUT_FASTA = "artificial_dna.fasta"

DNA_ALPHABET = "ACGT"
GENOME_LEN = random.randint(200, 400)   # 200-400 bp
NUM_TES = random.randint(3, 4)

# TE sequences (artificial, distinct)
TE_LIBRARY = {
    "TE1": "TACGTTGACCTGATGCTAGT",
    "TE2": "GGATCCGATGTTACCGGAA",
    "TE3": "CTAGGCTTACGATCGTAGC",
    "TE4": "AAGCTTGGGACCATTTGCA"
}

# optional: small terminal repeats (helps “TE feel”)
LEFT_REPEAT  = "TTA"
RIGHT_REPEAT = "TTA"


def rand_dna(n: int) -> str:
    return "".join(random.choice(DNA_ALPHABET) for _ in range(n))


def insert_non_overlapping(genome: str, inserts: list):
    """
    inserts: list of (name, seq) to insert
    returns: new_genome, placements list of dicts with start/end 1-based
    """
    placements = []
    occupied = []  # list of (start0, end0) occupied intervals in current genome coordinates after each insertion

    # We will insert sequentially; coordinates are handled by rebuilding the string.
    # To keep it simple, we choose an insertion position in the CURRENT genome each time,
    # ensuring no overlap with already inserted TE intervals.

    for te_name, te_seq in inserts:
        te_full = LEFT_REPEAT + te_seq + RIGHT_REPEAT
        te_len = len(te_full)

        # try multiple random positions
        for _ in range(5000):
            pos0 = random.randint(0, len(genome))  # insertion between bases
            # insertion interval becomes [pos0, pos0+te_len) in new genome
            # check overlap against already placed intervals (in current genome coords)
            ok = True
            for a, b in occupied:
                # if new interval overlaps old
                if not (pos0 + te_len <= a or pos0 >= b):
                    ok = False
                    break
            if ok:
                genome = genome[:pos0] + te_full + genome[pos0:]
                # update occupied: shift intervals after pos0 by te_len
                new_occupied = []
                for a, b in occupied:
                    if a >= pos0:
                        new_occupied.append((a + te_len, b + te_len))
                    else:
                        new_occupied.append((a, b))
                # add this TE interval
                new_occupied.append((pos0, pos0 + te_len))
                occupied = new_occupied

                placements.append({
                    "name": te_name,
                    "start": pos0 + 1,        # 1-based inclusive
                    "end": pos0 + te_len,     # 1-based inclusive
                    "length": te_len,
                    "sequence": te_full
                })
                break
        else:
            raise RuntimeError("Failed to place TE without overlap. Try again.")

    return genome, placements


def main():
    random.seed()  # random each run

    base_genome = rand_dna(GENOME_LEN)

    # pick 3-4 distinct TEs
    te_names = random.sample(list(TE_LIBRARY.keys()), k=NUM_TES)
    inserts = [(name, TE_LIBRARY[name]) for name in te_names]

    genome_with_tes, placements = insert_non_overlapping(base_genome, inserts)

    dataset = {
        "genome_length": len(genome_with_tes),
        "base_genome_length": GENOME_LEN,
        "tes_inserted": placements,
        "te_library": {
            name: LEFT_REPEAT + TE_LIBRARY[name] + RIGHT_REPEAT
            for name in te_names
        },
        "genome": genome_with_tes
    }

    # save JSON (for detection ground truth + TE library)
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(dataset, f, indent=2)

    # save FASTA
    with open(OUT_FASTA, "w", encoding="utf-8") as f:
        f.write(">artificial_dna_with_TEs\n")
        for i in range(0, len(genome_with_tes), 80):
            f.write(genome_with_tes[i:i+80] + "\n")

    print(f"[+] Generated genome length: {len(genome_with_tes)} bp")
    print(f"[+] Inserted {len(placements)} TEs")
    print(f"[+] Saved: {OUT_JSON}, {OUT_FASTA}")
    print("\nTE placements (ground truth):")
    for p in placements:
        print(f"  {p['name']}: {p['start']}-{p['end']} (len={p['length']})")


if __name__ == "__main__":
    main()
