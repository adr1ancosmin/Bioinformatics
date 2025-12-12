#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import math
import random
from typing import List, Dict

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw

# ---------------- basic setup ----------------
OUT_DIR = "lab6_ex2_output"
CUT_SITE = "GAATTC"           # EcoRI
N_GENOMES = 10
RNG = random.Random(77)

os.makedirs(OUT_DIR, exist_ok=True)

# -------------- DNA generation --------------
def make_dna(length: int, gc_bias: float, seed: int) -> str:
    rng = random.Random(seed)
    at = (1.0 - gc_bias) / 2.0
    gc = gc_bias / 2.0
    pool = (["A"] * int(at * 100) +
            ["T"] * int(at * 100) +
            ["G"] * int(gc * 100) +
            ["C"] * int(gc * 100))
    return "".join(rng.choice(pool) for _ in range(length))


def build_panel() -> Dict[str, str]:
    specs = [
        ("DNA01", 13520, 0.43),
        ("DNA02", 13280, 0.46),
        ("DNA03", 13840, 0.44),
        ("DNA04", 13410, 0.45),
        ("DNA05", 13660, 0.42),
        ("DNA06", 13190, 0.47),
        ("DNA07", 13750, 0.44),
        ("DNA08", 13330, 0.45),
        ("DNA09", 13570, 0.43),
        ("DNA10", 13240, 0.46),
    ]
    panel = {}
    for i, (name, length, gc) in enumerate(specs, start=1):
        panel[name] = make_dna(length, gc, seed=100 + i)
    return panel


# -------------- restriction digest --------------
def ecoRI_fragments(seq: str, site: str = CUT_SITE) -> List[int]:
    cuts = []
    pos = seq.find(site)
    while pos != -1:
        cuts.append(pos + 1)      # cut after the leading G
        pos = seq.find(site, pos + 1)
    positions = [0] + sorted(cuts) + [len(seq)]
    return [positions[i + 1] - positions[i] for i in range(len(positions) - 1)]


def ladder_sizes() -> List[int]:
    return [
        10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000,
        1500, 1200, 1000, 800, 700, 600, 500, 400,
        300, 250, 200, 150, 100,
    ]


# -------------- gel rendering --------------
def bp_to_y(bp: int, height: int,
           min_bp: int = 100, max_bp: int = 15000) -> int:
    bp = max(min_bp, min(max_bp, bp))
    log_min = math.log10(min_bp)
    log_max = math.log10(max_bp)
    t = (math.log10(bp) - log_min) / (log_max - log_min)
    t = 1.0 - t
    return int(t * (height - 1))


def draw_lane(fragments: List[int], width: int, height: int) -> Image.Image:
    img = Image.new("L", (width, height), 0)
    draw = ImageDraw.Draw(img)
    for bp in fragments:
        y = bp_to_y(bp, height)
        band_thick = max(1, int(1 + 5 * (1 - y / height)))
        draw.rectangle(
            [0, max(0, y - band_thick // 2), width - 1, min(height - 1, y + band_thick // 2)],
            fill=230,
        )
    return img


def make_gel(all_frags: Dict[str, List[int]]) -> str:
    labels = ["Ladder"] + list(all_frags.keys())
    lane_lists = [ladder_sizes()] + [all_frags[k] for k in all_frags]

    h = 900
    w_lane = 60
    spacer = 20
    text_margin = 150
    n_lanes = len(labels)

    total_w = n_lanes * w_lane + (n_lanes + 1) * spacer + text_margin
    gel = Image.new("L", (total_w, h), 0)
    draw = ImageDraw.Draw(gel)

    x = spacer
    for label, frags in zip(labels, lane_lists):
        lane = draw_lane(frags, w_lane, h)
        gel.paste(lane, (x, 0))

        lab_img = Image.new("L", (w_lane, 16), 0)
        ImageDraw.Draw(lab_img).text((2, 0), label, fill=255)
        lab_rot = lab_img.rotate(90, expand=1)
        gel.paste(lab_rot, (x + w_lane // 2 - lab_rot.width // 2, h - lab_rot.height - 4))

        x += w_lane + spacer

    # right-side size labels
    for bp in ladder_sizes():
        y = bp_to_y(bp, h)
        draw.text((total_w - text_margin + 30, y - 6), f"{bp} bp", fill=200)
        draw.line([(total_w - text_margin + 15, y), (total_w - text_margin + 25, y)], fill=180, width=1)

    path = os.path.join(OUT_DIR, "ecoRI_gel.png")
    gel.save(path)
    return path


# -------------- overlay plot --------------
def make_overlay_plot(all_frags: Dict[str, List[int]], richest_name: str, richest_count: int) -> str:
    plt.figure(figsize=(10, 6))
    for name, frags in all_frags.items():
        sizes = sorted(frags, reverse=True)
        xs = np.arange(1, len(sizes) + 1)
        plt.plot(xs, sizes, marker="o", linewidth=1.4, label=name)

    plt.xlabel("Fragment index (largest → smallest)")
    plt.ylabel("Fragment size (bp)")
    plt.title(f"EcoRI fragment sizes — most bands: {richest_name} ({richest_count})")
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()

    path = os.path.join(OUT_DIR, "ecoRI_fragments_overlay.png")
    plt.savefig(path, dpi=200)
    plt.close()
    return path


# -------------- main --------------
def main():
    genomes = build_panel()

    digest_map: Dict[str, List[int]] = {}
    for name, seq in genomes.items():
        frags = ecoRI_fragments(seq)
        digest_map[name] = frags

    counts = {name: len(flist) for name, flist in digest_map.items()}
    richest_name = max(counts, key=counts.get)
    richest_count = counts[richest_name]

    print("EcoRI digest summary:")
    for name in sorted(digest_map):
        total_len = len(genomes[name])
        bands = counts[name]
        print(f"  {name}: {bands} fragments, total length {total_len} bp")

    print(f"\nMost bands: {richest_name} ({richest_count} fragments)")

    gel_path = make_gel(digest_map)
    overlay_path = make_overlay_plot(digest_map, richest_name, richest_count)

    summary_path = os.path.join(OUT_DIR, "ecoRI_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write("EcoRI digest results\n")
        fh.write("====================\n\n")
        for name in sorted(digest_map):
            sizes = sorted(digest_map[name], reverse=True)
            fh.write(f"{name}  (len {len(genomes[name])} bp)\n")
            fh.write(f"  bands: {len(sizes)}\n")
            fh.write(f"  sizes (bp): {', '.join(str(s) for s in sizes)}\n\n")
        fh.write(f"Most fragments: {richest_name} ({richest_count})\n")
        fh.write(f"Gel image: {gel_path}\n")
        fh.write(f"Overlay plot: {overlay_path}\n")

    print("\nFiles written in:", OUT_DIR)
    print("  -", os.path.basename(gel_path))
    print("  -", os.path.basename(overlay_path))
    print("  -", os.path.basename(summary_path))


if __name__ == "__main__":
    main()
