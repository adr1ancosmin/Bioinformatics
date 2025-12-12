# restriction digest + fake gel image
# needs: pip install matplotlib

import matplotlib.pyplot as plt

enzymes = {
    "EcoRI":  ("GAATTC", 1),  # G|AATTC
    "BamHI":  ("GGATCC", 1),  # G|GATCC
    "HindIII":("AAGCTT", 1),  # A|AGCTT
    "TaqI":   ("TCGA", 1),    # T|CGA
    "HaeIII": ("GGCC",  2),   # GG|CC
}

# simple size marker ladder in bp
ladder_sizes = [1000, 800, 600, 400, 300, 200, 100]

def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line.upper())
    return "".join(seq)

def find_cuts(seq, motif, cut_pos):
    pos = []
    i = 0
    while True:
        i = seq.find(motif, i)
        if i == -1:
            break
        pos.append(i + cut_pos)   # cut inside motif
        i += 1                    # move one base
    return sorted(pos)

def get_fragments(seq_len, cuts):
    if not cuts:
        return [seq_len]
    cuts = sorted(cuts)
    pts = [0] + cuts + [seq_len]
    frags = []
    for i in range(len(pts) - 1):
        frags.append(pts[i+1] - pts[i])
    return frags

def draw_gel_image(frag_by_enzyme, ladder, filename="gel.png"):
    lanes = list(frag_by_enzyme.keys()) + ["Marker"]

    # collect all sizes for scaling
    all_sizes = []
    for frags in frag_by_enzyme.values():
        all_sizes.extend(frags)
    all_sizes.extend(ladder)
    max_size = max(all_sizes)
    min_size = min(all_sizes)

    # map size to y position (big = near top, small = bottom)
    def size_to_y(size):
        frac = (size - min_size) / (max_size - min_size + 1e-9)
        return 1.0 - 0.9 * frac   # keep a bit of margin at top/bottom

    fig, ax = plt.subplots(figsize=(4, 6))
    ax.set_facecolor("black")
    fig.patch.set_facecolor("black")

    # draw one rectangular well row
    for i in range(len(lanes)):
        ax.add_patch(plt.Rectangle((i + 0.1, 0.95), 0.8, 0.01, color="grey"))

    # draw bands for each lane
    for lane_idx, name in enumerate(lanes):
        if name == "Marker":
            sizes = ladder
            color = "white"
        else:
            sizes = frag_by_enzyme[name]
            color = "white"

        for s in sizes:
            y = size_to_y(s)
            # band thickness depends a bit on size (optional)
            band_height = 0.01
            ax.add_patch(
                plt.Rectangle(
                    (lane_idx + 0.15, y - band_height / 2),
                    0.7,
                    band_height,
                    color=color,
                )
            )

    # formatting
    ax.set_xlim(0, len(lanes))
    ax.set_ylim(0, 1)
    ax.invert_yaxis()  # top = wells
    ax.set_xticks([i + 0.5 for i in range(len(lanes))])
    ax.set_xticklabels([name for name in lanes], color="white")
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"gel image saved as {filename}")

def main():
    seq = read_fasta("sequence.fasta")
    L = len(seq)
    print(f"sequence length: {L} bp\n")

    frag_by_enzyme = {}

    for name, (motif, cut) in enzymes.items():
        cuts = find_cuts(seq, motif, cut)
        frags = get_fragments(L, cuts)
        frag_by_enzyme[name] = frags

        print(f"=== {name} ({motif}) ===")
        print("number of cuts:", len(cuts))
        print("fragment lengths (bp):", frags)
        print()

    draw_gel_image(frag_by_enzyme, ladder_sizes, filename="gel.png")

if __name__ == "__main__":
    main()
