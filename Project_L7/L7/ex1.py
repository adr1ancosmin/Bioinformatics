#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
from pathlib import Path
from typing import Dict, List, Tuple
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np

# ----------------------------- CONFIG -----------------------------
OUTDIR = "lab7_ex1_output"
N_GENOMES = 11                        # 1 arbitrary + 10 “influenza-like”
K_MIN = 2
K_MAX = 6
MIN_REPEAT_BLOCK = 2
SEQ_MIN = 1500
SEQ_MAX = 2900
SEED = 101
TOP_UNITS = 30                         # top-N repeat units in bar plot
random.seed(SEED)

# ---------------------- DNA GENERATION ----------------------
def make_dna(length: int, gc_bias: float, seed: int) -> str:
    rng = random.Random(seed)
    at = (1.0 - gc_bias) / 2
    gc = gc_bias / 2
    pool = (
        ["A"] * int(at * 100)
        + ["T"] * int(at * 100)
        + ["G"] * int(gc * 100)
        + ["C"] * int(gc * 100)
    )
    return "".join(rng.choice(pool) for _ in range(length))


def build_dna_panel() -> Dict[str, str]:
    genomes = {}
    genomes["custom_01"] = make_dna(
        random.randint(SEQ_MIN, SEQ_MAX), 0.43, seed=500
    )

    for i in range(1, 11):
        length = random.randint(SEQ_MIN, SEQ_MAX)
        gc = 0.40 + random.uniform(-0.03, 0.03)
        genomes[f"virus_{i:02d}"] = make_dna(length, gc, seed=600 + i)

    return genomes


# ---------------------- TANDEM REPEAT DETECTION ----------------------
def count_tandem_blocks(seq: str, kmin: int, kmax: int, min_block: int) -> Counter:
    seq = seq.upper()
    n = len(seq)
    counts = Counter()
    i = 0

    while i < n:
        moved = False
        for k in range(kmax, kmin - 1, -1):
            if i + k * min_block > n:
                continue

            unit = seq[i : i + k]
            m = 1
            while i + (m + 1) * k <= n and seq[i + m * k : i + (m + 1) * k] == unit:
                m += 1

            if m >= min_block:
                counts[unit] += 1
                i += m * k
                moved = True
                break

        if not moved:
            i += 1

    return counts


# ---------------------- PLOTTING ----------------------
def plot_grouped(series: List[Tuple[str, Counter]], outfile: Path):
    total = Counter()
    for _, cnt in series:
        total.update(cnt)

    units = [u for u, _ in total.most_common(TOP_UNITS)]
    if not units:
        plt.figure(figsize=(10, 4))
        plt.title("No tandem repeats detected")
        plt.savefig(outfile, dpi=220)
        plt.close()
        return

    x = np.arange(len(units))
    width = 0.75 / len(series)

    plt.figure(figsize=(max(12, 0.6 * len(units)), 6))
    for i, (label, cnt) in enumerate(series):
        values = [cnt.get(u, 0) for u in units]
        plt.bar(x + i * width, values, width=width, label=label)

    plt.xticks(x + (len(series) - 1) * width / 2, units, rotation=90)
    plt.xlabel("repeat unit (tandem)")
    plt.ylabel("frequency (blocks)")
    plt.title(f"Tandem repeats across custom genomes (k={K_MIN}..{K_MAX}, min×{MIN_REPEAT_BLOCK})")
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(outfile, dpi=220)
    plt.close()


# ---------------------- MAIN ----------------------
def main():
    out_dir = Path(OUTDIR)
    out_dir.mkdir(exist_ok=True)

    genomes = build_dna_panel()

    series: List[Tuple[str, Counter]] = []
    for name, seq in genomes.items():
        repeat_profile = count_tandem_blocks(
            seq, kmin=K_MIN, kmax=K_MAX, min_block=MIN_REPEAT_BLOCK
        )
        series.append((name, repeat_profile))

    out_path = out_dir / "tandem_repeats_custom.png"
    plot_grouped(series, out_path)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
