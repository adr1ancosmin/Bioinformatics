#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import json
import random
import textwrap
from collections import defaultdict, Counter

import numpy as np
import matplotlib.pyplot as plt

N_GENOMES = 10
READ_LEN = 150
COVERAGE = 20
K = 31
RANDOM_SEED = 42

random.seed(RANDOM_SEED)


def build_sequence(length, motif):
    reps = length // len(motif) + 1
    return (motif * reps)[:length]


def create_dna_panel():
    configs = [
        ("V01", 2610, "ATGCGT"),
        ("V02", 2560, "GCCGGC"),
        ("V03", 2565, "ATCCGC"),
        ("V04", 5852, "GCGGAT"),
        ("V05", 4380, "ATATGC"),
        ("V06", 1930, "ATGCGC"),
        ("V07", 1935, "ATATGC"),
        ("V08", 2079, "GCGTGC"),
        ("V09", 1680, "ATATAT"),
        ("V10", 2734, "GCCGAT"),
    ]
    genomes = []
    for name, length, motif in configs:
        seq = build_sequence(length, motif)
        header = f"{name}|DNA"
        genomes.append({"label": name, "header": header, "seq": seq})
    return genomes


def gc_percent(seq):
    gc = sum(1 for b in seq if b in ("G", "C"))
    at = sum(1 for b in seq if b in ("A", "T"))
    tot = gc + at
    return 100.0 * gc / tot if tot else 0.0


def make_reads(seq, read_len, coverage):
    target = max(1, int(round(coverage * len(seq) / read_len)))
    reads = []
    if len(seq) <= read_len:
        return [seq] * target
    for _ in range(target):
        start = random.randint(0, len(seq) - read_len)
        reads.append(seq[start:start + read_len])
    return reads


def make_graph(reads, k):
    edges = defaultdict(Counter)
    indeg = Counter()
    outdeg = Counter()
    for r in reads:
        if len(r) < k:
            continue
        for i in range(len(r) - k + 1):
            frag = r[i:i + k]
            if "N" in frag:
                continue
            a = frag[:-1]
            b = frag[1:]
            edges[a][b] += 1
            outdeg[a] += 1
            indeg[b] += 1
    return edges, indeg, outdeg


def walk_contigs(graph):
    contigs = []
    visited = set()
    for start in list(graph.keys()):
        if start in visited:
            continue
        cur = start
        seq = cur
        visited.add(cur)
        while graph.get(cur):
            nxt = graph[cur].most_common(1)[0][0]
            seq += nxt[-1]
            cur = nxt
            if cur in visited:
                break
            visited.add(cur)
            if len(seq) > 2_000_000:
                break
        contigs.append(seq)
    return contigs


def run_assembly(reads, k):
    t0 = time.perf_counter()
    graph, indeg, outdeg = make_graph(reads, k)
    contigs = walk_contigs(graph)
    elapsed_ms = (time.perf_counter() - t0) * 1000.0
    lengths = sorted((len(c) for c in contigs), reverse=True)
    n50 = 0
    if lengths:
        half = sum(lengths) / 2.0
        acc = 0
        for L in lengths:
            acc += L
            if acc >= half:
                n50 = L
                break
    return elapsed_ms, n50, len(contigs)


def format_bp(n):
    return f"{n:,} bp"


def main():
    os.makedirs("out/genomes", exist_ok=True)
    os.makedirs("out/reads", exist_ok=True)
    os.makedirs("out/results", exist_ok=True)

    panel = create_dna_panel()
    results = []

    for idx, item in enumerate(panel, start=1):
        label = f"V{idx:02d}"
        header = item["header"]
        seq = item["seq"]

        fasta_path = f"out/genomes/{label}.fasta"
        with open(fasta_path, "w") as f:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")

        gc = gc_percent(seq)
        reads = make_reads(seq, READ_LEN, COVERAGE)

        fastq_path = f"out/reads/{label}_reads.fq"
        with open(fastq_path, "w") as f:
            for j, r in enumerate(reads, start=1):
                f.write(f"@{label}_r{j}\n{r}\n+\n{'I' * len(r)}\n")

        ms, n50, n_contigs = run_assembly(reads, K)

        rec = {
            "index": idx,
            "label": label,
            "header": header,
            "genome_len": len(seq),
            "gc_percent": gc,
            "assembly_time_ms": ms,
            "n50_approx": n50,
            "n_contigs": n_contigs,
            "n_reads": len(reads),
        }
        results.append(rec)
        print(
            f"[{label}] len={len(seq)} bp | GC={gc:.2f} % | reads={len(reads)} | "
            f"time={ms:.1f} ms | N50≈{n50} | contigs={n_contigs}"
        )

    summary_path = "out/results/summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    xs = [r["gc_percent"] for r in results]
    ys = [r["assembly_time_ms"] for r in results]
    labels = [r["label"] for r in results]
    lens_bp = [r["genome_len"] for r in results]

    plt.figure(figsize=(9.6, 6.0))
    plt.scatter(xs, ys)
    for x, y, label, L in zip(xs, ys, labels, lens_bp):
        plt.text(x, y, f"{label} • {format_bp(L)}", fontsize=8, ha="left", va="bottom")
    plt.xlabel("Overall GC (%)")
    plt.ylabel("Assembly time (ms)")
    plt.title("Toy assembly timing vs GC% (10 DNA sequences)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    gc_arr = np.array(xs, dtype=float)
    time_arr = np.array(ys, dtype=float)
    corr_gc_time = float(np.corrcoef(gc_arr, time_arr)[0, 1]) if len(xs) > 1 else float("nan")

    len_arr = np.array([r["genome_len"] for r in results], dtype=float)
    corr_len_time = float(np.corrcoef(len_arr, time_arr)[0, 1]) if len(xs) > 1 else float("nan")

    fastest = min(results, key=lambda r: r["assembly_time_ms"])
    slowest = max(results, key=lambda r: r["assembly_time_ms"])
    leftmost = min(results, key=lambda r: r["gc_percent"])
    rightmost = max(results, key=lambda r: r["gc_percent"])

    lines = []
    lines.append("Point Distribution on the GC% vs Time Plot")
    lines.append("=" * 44)
    lines.append("")
    lines.append(f"Number of sequences: {len(results)}")
    lines.append(f"GC% range: {min(xs):.2f}% – {max(xs):.2f}%")
    lines.append(f"Assembly time range: {min(ys):.1f} ms – {max(ys):.1f} ms")
    lines.append(f"Correlation (GC% vs time): {corr_gc_time:.3f}")
    lines.append(f"Correlation (length vs time): {corr_len_time:.3f}")
    lines.append("")
    lines.append("Examples of individual points")
    lines.append("-" * 31)
    lines.append(
        f"- Fastest point: {fastest['label']} "
        f"(time {fastest['assembly_time_ms']:.1f} ms, GC {fastest['gc_percent']:.2f}%, "
        f"length {fastest['genome_len']} bp)."
    )
    lines.append(
        f"- Slowest point: {slowest['label']} "
        f"(time {slowest['assembly_time_ms']:.1f} ms, GC {slowest['gc_percent']:.2f}%, "
        f"length {slowest['genome_len']} bp)."
    )
    lines.append(
        f"- Lowest GC (left side of X-axis): {leftmost['label']} "
        f"(GC {leftmost['gc_percent']:.2f}%, time {leftmost['assembly_time_ms']:.1f} ms)."
    )
    lines.append(
        f"- Highest GC (right side of X-axis): {rightmost['label']} "
        f"(GC {rightmost['gc_percent']:.2f}%, time {rightmost['assembly_time_ms']:.1f} ms)."
    )
    lines.append("")
    lines.append("How to read the plot")
    lines.append("-" * 20)
    lines.append(
        textwrap.fill(
            "The horizontal axis shows GC%. Sequences with lower GC content are placed more to the left, "
            "and those with higher GC content appear to the right. The vertical axis shows the time needed "
            "for the small de Bruijn assembly run on reads sampled from each sequence.",
            width=90,
        )
    )
    lines.append("")
    lines.append(
        textwrap.fill(
            "Points that share a similar GC% but sit at different heights usually differ in length or graph "
            "complexity: a longer or more repetitive sequence generates more k-mers and tends to require more "
            "work to traverse, so the point moves upward. Points with similar assembly time but different GC% "
            "indicate that GC composition alone is not the main factor controlling runtime in this toy setup.",
            width=90,
        )
    )
    lines.append("")
    lines.append("Summary table")
    lines.append("-" * 13)
    lines.append(
        "Label  GC%    Time(ms)  Length(bp)  Reads  ~N50  Contigs  Header"
    )
    lines.append(
        "-----  -----  --------  ----------  -----  ----  -------  ------"
    )
    for r in results:
        lines.append(
            f"{r['label']:>4}  {r['gc_percent']:5.2f}  "
            f"{r['assembly_time_ms']:8.1f}  {r['genome_len']:10d}  "
            f"{r['n_reads']:5d}  {r['n50_approx']:4d}  {r['n_contigs']:7d}  {r['header']}"
        )

    report_text = "\n".join(lines)
    print("\n" + report_text + "\n")

    report_path = "out/results/explanation.txt"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_text)

    print("Done.")
    print(f"- Summary JSON: {summary_path}")
    print(f"- Text report:  {report_path}")
    print("Outputs are in ./out/")

if __name__ == "__main__":
    main()
