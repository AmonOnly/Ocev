#!/usr/bin/env python3
"""
Simple script to read evolution_stats.csv and plot max and average fitness over generations.
Saves output as evolution_plot.png in workspace root.
"""
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CSV_GLOB = ROOT.glob('evolution_stats*.csv')
OUT_BEST = ROOT / 'evolution_best.png'
OUT_STD = ROOT / 'evolution_std.png'


def main():
    # find all matching CSVs (support single file or multiple runs)
    files = sorted([p for p in ROOT.glob('evolution_stats*.csv')])
    if not files:
        print(f"No evolution_stats CSV files found in {ROOT}")
        return 1

    # read all dataframes into a list
    dfs = []
    for p in files:
        df = pd.read_csv(p)
        if 'generation' not in df.columns:
            df = df.reset_index().rename(columns={'index': 'generation'})
        dfs.append(df)

    # determine maximum generation across runs
    max_gen = max(df['generation'].max() for df in dfs)

    # align all runs to same generation index (0..max_gen)
    aligned = {}
    for i, df in enumerate(dfs):
        s_max = pd.Series(index=range(0, int(max_gen) + 1), dtype=float)
        s_avg = pd.Series(index=range(0, int(max_gen) + 1), dtype=float)
        s_max[:] = float('nan')
        s_avg[:] = float('nan')
        for _, row in df.iterrows():
            g = int(row['generation'])
            if 'max_fitness' in row:
                s_max[g] = row['max_fitness']
            if 'avg_fitness' in row:
                s_avg[g] = row['avg_fitness']
        aligned[i] = {'max': s_max, 'avg': s_avg}

    # build DataFrames for aggregation
    max_df = pd.DataFrame({i: aligned[i]['max'] for i in aligned})
    avg_df = pd.DataFrame({i: aligned[i]['avg'] for i in aligned})

    # compute per-generation statistics across runs
    max_best = max_df.max(axis=1)
    max_std = max_df.std(axis=1)
    avg_best = avg_df.max(axis=1)
    avg_std = avg_df.std(axis=1)

    # prefer seaborn style if available, otherwise fallback to a built-in style
    try:
        plt.style.use('seaborn-darkgrid')
    except Exception:
        try:
            plt.style.use('seaborn')
        except Exception:
            plt.style.use('classic')
    gens = max_best.index.astype(int)

    # Plot 1: best (max across runs) for max_fitness and avg_fitness
    fig1, ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(gens, max_best, label='Max Fitness (best across runs)', color='tab:blue', marker='o')
    if not avg_best.isna().all():
        ax1.plot(gens, avg_best, label='Avg Fitness (best across runs)', color='tab:orange', marker='o')
    ax1.set_xlabel('Generation')
    ax1.set_ylabel('Fitness (number of satisfied clauses)')
    ax1.set_title('Best Fitness per Generation (across runs)')
    ax1.legend()
    fig1.tight_layout()
    fig1.savefig(OUT_BEST, dpi=200)
    print(f"Saved best plot to {OUT_BEST}")

    # Plot 2: standard deviation across runs for max_fitness and avg_fitness
    fig2, ax2 = plt.subplots(figsize=(10, 5))
    ax2.plot(gens, max_std, label='StdDev of Max Fitness', color='tab:green', marker='o')
    if not avg_std.isna().all():
        ax2.plot(gens, avg_std, label='StdDev of Avg Fitness', color='tab:red', marker='o')
    ax2.set_xlabel('Generation')
    ax2.set_ylabel('Standard Deviation')
    ax2.set_title('Standard Deviation of Fitness per Generation (across runs)')
    ax2.legend()
    fig2.tight_layout()
    fig2.savefig(OUT_STD, dpi=200)
    print(f"Saved stddev plot to {OUT_STD}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
