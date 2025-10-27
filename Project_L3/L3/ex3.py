import os
import math
from collections import Counter
from typing import Iterator, Tuple, List, Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

FASTA_PATH = "dna.fasta"  
WINDOW     = 9              
STEP       = 1              
OUTDIR     = "results"      
PREFIX     = None           
NA_MOLAR   = 0.001          

THRESHOLDS: Dict[str, float] = {
    "Tm_simple":  25.0,
    "Tm_complex": -10.0,
}

def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    header = None
    seq_chunks: List[str] = []
    with open(path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_chunks).upper()
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header, ''.join(seq_chunks).upper()

def computeMeltingTempSimple(seq: str) -> float:
    s = ''.join(ch for ch in seq.upper() if ch in 'ACGT')
    if not s:
        return float('nan')
    A = s.count("A"); C = s.count("C"); G = s.count("G"); T = s.count("T")
    return 4 * (C + G) + 2 * (A + T)

def computeMeltingTempComplex(seq: str, na_molar: float = NA_MOLAR) -> float:
    s = ''.join(ch for ch in seq.upper() if ch in 'ACGT')
    if not s:
        return float('nan')
    if na_molar <= 0:
        raise ValueError("[Na+] must be positive (mol/L).")
    length = len(s)
    gc_pct = 100.0 * (s.count("G") + s.count("C")) / length
    return 81.5 + 16.6 * math.log10(na_molar) + 0.41 * gc_pct - (600.0 / length)

def window_metrics(seq: str, window: int, step: int = 1, na_molar: float = NA_MOLAR) -> pd.DataFrame:
    n = len(seq)
    data = {
        'win': [],
        'position': [],
        'A': [], 'C': [], 'G': [], 'T': [],
        'Tm_simple': [], 'Tm_complex': []
    }
    valid = set('ACGT')

    win_id = 0
    for start in range(0, n - window + 1, step):
        win_id += 1
        wseq = seq[start:start + window]
        filtered = ''.join(b for b in wseq if b in valid)
        counts = Counter(filtered)
        denom = sum(counts.values())
        if denom == 0:
            a = c = g = t = math.nan
            tm_s = tm_c = math.nan
        else:
            a = 100.0 * counts.get('A', 0) / denom
            c = 100.0 * counts.get('C', 0) / denom
            g = 100.0 * counts.get('G', 0) / denom
            t = 100.0 * counts.get('T', 0) / denom
            tm_s = computeMeltingTempSimple(filtered)
            tm_c = computeMeltingTempComplex(filtered, na_molar)
        center = start + (window // 2) + 1

        data['win'].append(win_id)
        data['position'].append(center)
        data['A'].append(a); data['C'].append(c); data['G'].append(g); data['T'].append(t)
        data['Tm_simple'].append(tm_s); data['Tm_complex'].append(tm_c)

    return pd.DataFrame(data)

def area_above_threshold(x: np.ndarray, y: np.ndarray, thr: float) -> float:
    above = np.clip(y - thr, 0, None)
    return float(np.trapz(above, x))

def plot_tm_signals_with_thresholds(df: pd.DataFrame,
                                    thresholds: Dict[str, float],
                                    title: str,
                                    out_png: str) -> Dict[str, Tuple[float, int]]:
    x = df['position'].to_numpy(dtype=float)
    y_simple  = df['Tm_simple'].to_numpy(dtype=float)
    y_complex = df['Tm_complex'].to_numpy(dtype=float)

    thr_simple  = thresholds['Tm_simple']
    thr_complex = thresholds['Tm_complex']

    plt.figure(figsize=(10, 5), dpi=120)
    plt.plot(x, y_simple,  label="Tm_simple",  linewidth=2.0, color='tab:blue')
    plt.plot(x, y_complex, label="Tm_complex", linewidth=2.0, color='tab:orange')

    plt.axhline(thr_simple,  linestyle='--', linewidth=1.5, color='tab:blue',   label=f"threshold simple = {thr_simple:g} °C")
    plt.axhline(thr_complex, linestyle='--', linewidth=1.5, color='tab:orange', label=f"threshold complex = {thr_complex:g} °C")

    stats = {}
    signals = {'Tm_simple': (y_simple, 'tab:blue'), 'Tm_complex': (y_complex, 'tab:orange')}
    
    for signal_name, (y_data, color) in signals.items():
        valid_mask = ~np.isnan(y_data)
        if np.any(valid_mask):
            valid_y = y_data[valid_mask]
            valid_x = x[valid_mask]
            valid_indices = np.where(valid_mask)[0]
            
            min_idx = np.nanargmin(valid_y)
            min_val = valid_y[min_idx]
            min_pos = valid_x[min_idx]
            min_win = df.iloc[valid_indices[min_idx]]['win']
            
            max_idx = np.nanargmax(valid_y)
            max_val = valid_y[max_idx]
            max_pos = valid_x[max_idx]
            max_win = df.iloc[valid_indices[max_idx]]['win']
            
            plt.scatter([min_pos], [min_val], s=50, color=color, marker='o', zorder=5)
            plt.scatter([max_pos], [max_val], s=50, color=color, marker='s', zorder=5)
            
            plt.annotate(f'min {signal_name}={min_val:.1f}', 
                        xy=(min_pos, min_val), xytext=(5, 5), 
                        textcoords='offset points', fontsize=8, 
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
            plt.annotate(f'max {signal_name}={max_val:.1f}', 
                        xy=(max_pos, max_val), xytext=(5, -15), 
                        textcoords='offset points', fontsize=8,
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
            
            stats[f'{signal_name}_min'] = (float(min_val), int(min_win))
            stats[f'{signal_name}_max'] = (float(max_val), int(max_win))

    plt.xlabel("Window center position (bp)")
    plt.ylabel("Melting temperature (°C)")
    plt.title(title)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(out_png)
    plt.show()
    plt.close()
    
    return stats

def _true_runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    runs = []
    start = None
    N = len(mask)
    for i in range(N):
        if mask[i] and start is None:
            start = i + 1
        elif (not mask[i]) and start is not None:
            runs.append((start, i))
            start = None
    if start is not None:
        runs.append((start, N))
    return runs

def plot_regions_above_threshold(df: pd.DataFrame,
                                 thresholds: Dict[str, float],
                                 out_png: str) -> None:
    wins = df['win'].to_numpy()
    y_simple  = df['Tm_simple'].to_numpy(dtype=float)
    y_complex = df['Tm_complex'].to_numpy(dtype=float)

    thr_simple  = thresholds['Tm_simple']
    thr_complex = thresholds['Tm_complex']

    mask_s = y_simple  > thr_simple
    mask_c = y_complex > thr_complex

    runs_s = _true_runs(mask_s)
    runs_c = _true_runs(mask_c)

    ranges_s = [(start - 0.5, (end - start + 1)) for start, end in runs_s]
    ranges_c = [(start - 0.5, (end - start + 1)) for start, end in runs_c]

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 4), dpi=120)

    ax1.broken_barh(ranges_s, (0.25, 0.5), facecolors='tab:blue')
    ax1.set_ylim(0, 1)
    ax1.set_yticks([0.5])
    ax1.set_yticklabels(['P1'])
    ax1.set_title(f"P1 - Regions Above Threshold ({thr_simple:.1f}°C)")
    ax1.grid(False)

    ax2.broken_barh(ranges_c, (0.25, 0.5), facecolors='tab:orange')
    ax2.set_ylim(0, 1)
    ax2.set_yticks([0.5])
    ax2.set_yticklabels(['P2'])
    ax2.set_title(f"P2 - Regions Above Threshold ({thr_complex:.1f}°C)")
    ax2.set_xlabel("Window Number")
    ax2.grid(False)

    nwin = int(wins.max()) if len(wins) else 0
    ax2.set_xlim(0, max(1, nwin) + 1)

    plt.tight_layout()
    plt.savefig(out_png)
    plt.show()
    plt.close()

def sanitize_name(name: str) -> str:
    safe = ''.join(ch if ch.isalnum() or ch in ('-', '_') else '_' for ch in name)
    return safe[:80] or 'sequence'

def main():
    if WINDOW <= 0:
        raise SystemExit("WINDOW must be positive")
    if STEP <= 0:
        raise SystemExit("STEP must be positive")
    if NA_MOLAR <= 0:
        raise SystemExit("NA_MOLAR must be positive (mol/L)")
    for key in ("Tm_simple", "Tm_complex"):
        if key not in THRESHOLDS:
            raise SystemExit(f"THRESHOLDS must define '{key}'")
    if not os.path.exists(FASTA_PATH):
        raise SystemExit(f"FASTA not found: {FASTA_PATH}")

    os.makedirs(OUTDIR, exist_ok=True)
    summary_rows = []

    for idx, (hdr, seq) in enumerate(read_fasta(FASTA_PATH), start=1):
        if len(seq) < WINDOW:
            print(f"[skip] Sequence {idx} ('{hdr}') shorter than window ({len(seq)} < {WINDOW})")
            continue

        df = window_metrics(seq, WINDOW, STEP, NA_MOLAR)
        base_prefix = PREFIX or sanitize_name(hdr or f"seq{idx}")

        chart1_png = os.path.join(OUTDIR, f"{base_prefix}_signals_with_thresholds.png")
        stats = plot_tm_signals_with_thresholds(
            df, THRESHOLDS,
            title=f"Tm signals with thresholds (W={WINDOW}, step={STEP}, [Na+]={NA_MOLAR} M)\n{hdr}",
            out_png=chart1_png
        )

        chart2_png = os.path.join(OUTDIR, f"{base_prefix}_regions_above_threshold.png")
        plot_regions_above_threshold(df, THRESHOLDS, chart2_png)

        print(f"[ok] {hdr} | Tm_simple min={stats['Tm_simple_min'][0]:.2f} (win {stats['Tm_simple_min'][1]}), "
              f"max={stats['Tm_simple_max'][0]:.2f} (win {stats['Tm_simple_max'][1]}); "
              f"Tm_complex min={stats['Tm_complex_min'][0]:.2f} (win {stats['Tm_complex_min'][1]}), "
              f"max={stats['Tm_complex_max'][0]:.2f} (win {stats['Tm_complex_max'][1]})")

        print(f"[files] {chart1_png}\n        {chart2_png}")

        x = df['position'].to_numpy(dtype=float)
        a_simple  = area_above_threshold(x, df['Tm_simple'].to_numpy(dtype=float), THRESHOLDS["Tm_simple"])
        a_complex = area_above_threshold(x, df['Tm_complex'].to_numpy(dtype=float), THRESHOLDS["Tm_complex"])
        summary_rows.append({
            "sequence": hdr,
            "window": WINDOW,
            "step": STEP,
            "na_molar": NA_MOLAR,
            "threshold_simple": THRESHOLDS["Tm_simple"],
            "threshold_complex": THRESHOLDS["Tm_complex"],
            "area_above_simple": a_simple,
            "area_above_complex": a_complex,
            "Tm_simple_min": stats['Tm_simple_min'][0],
            "Tm_simple_min_win": stats['Tm_simple_min'][1],
            "Tm_simple_max": stats['Tm_simple_max'][0],
            "Tm_simple_max_win": stats['Tm_simple_max'][1],
            "Tm_complex_min": stats['Tm_complex_min'][0],
            "Tm_complex_min_win": stats['Tm_complex_min'][1],
            "Tm_complex_max": stats['Tm_complex_max'][0],
            "Tm_complex_max_win": stats['Tm_complex_max'][1],
        })

if __name__ == "__main__":
    main()
