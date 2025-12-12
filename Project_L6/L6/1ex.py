#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import random
import matplotlib.pyplot as plt

# ----------------- configuration -----------------
N_FRAGMENTS = 10
MIN_LEN = 100
MAX_LEN = 3000
MIN_TOTAL = 1000
MAX_TOTAL = 3000
RANDOM_SEED = 7
# -------------------------------------------------


def get_dna():
    """Return a header + one DNA sequence."""
    header = "DNA_01"
    motif = "ATGCGT"
    target_len = 2083
    reps = target_len // len(motif) + 1
    seq = (motif * reps)[:target_len]
    return header, seq


def check_length(seq, low, high):
    n = len(seq)
    if n < low:
        raise SystemExit(f"Sequence is too short: {n} nt (min {low}).")
    if n > high:
        raise SystemExit(f"Sequence is too long: {n} nt (max {high}).")
    print(f"Sequence length OK: {n} nt (within {low}-{high}).")
    return n


def sample_fragments(seq, n_frag, min_len, max_len):
    L = len(seq)
    max_len = min(max_len, L)
    if max_len < min_len:
        raise SystemExit("Fragment length range is incompatible with sequence length.")
    out = []
    for _ in range(n_frag):
        size = random.randint(min_len, max_len)
        start = random.randint(0, L - size)
        end = start + size
        out.append((start, end, size))
    return out


def lengths_to_positions(lengths):
    logs = [math.log10(x) for x in lengths]
    lo, hi = min(logs), max(logs)
    if hi == lo:
        return [0.5] * len(lengths)
    top, bottom = 0.1, 0.9
    span = bottom - top
    return [top + (hi - lg) / (hi - lo) * span for lg in logs]


def draw_gel(lengths, lane_label, title):
    positions = lengths_to_positions(lengths)

    fig, ax = plt.subplots(figsize=(4.5, 6.5))
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # outer frame
    ax.plot([0.18, 0.82, 0.82, 0.18, 0.18],
            [0.05, 0.05, 0.95, 0.95, 0.05],
            color="white", linewidth=1.6)

    # inner lane
    ax.plot([0.40, 0.60, 0.60, 0.40, 0.40],
            [0.07, 0.07, 0.93, 0.93, 0.07],
            color="white", linewidth=1.0)

    for y, L in sorted(zip(positions, lengths), key=lambda t: t[0]):
        ax.plot([0.41, 0.59], [y, y], linewidth=4.0, color="yellow")
        ax.text(0.62, y, f"{L} bp", va="center",
                ha="left", fontsize=9, color="white")

    ax.text(0.50, 0.03, lane_label, ha="center",
            va="center", fontsize=10, color="white")
    ax.set_title(title, fontsize=11, color="white", pad=10)

    plt.tight_layout()
    plt.show()


def main():
    random.seed(RANDOM_SEED)

    header, seq = get_dna()
    print(f"Sequence: {header}")
    check_length(seq, MIN_TOTAL, MAX_TOTAL)

    samples = sample_fragments(seq, N_FRAGMENTS, MIN_LEN, MAX_LEN)

    print("\nRandom samples (start..end, length):")
    for i, (s, e, size) in enumerate(samples, start=1):
        print(f"  {i:2d}: {s}..{e} (len={size} bp)")

    lengths = [size for (_s, _e, size) in samples]
    print("\nRandom sample lengths (bp):", lengths)

    draw_gel(
        lengths,
        lane_label=f"{N_FRAGMENTS} random samples",
        title=f"Fragments from {header}",
    )


if __name__ == "__main__":
    main()
