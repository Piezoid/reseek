#!/usr/bin/env python3

import sys
import os
import numpy as np
import bca_reader
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter
from collections import Counter

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 compute_product_alphabet.py <filename.bca>")
        sys.exit(1)
    
    filename = sys.argv[1]
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)
    
    bca = bca_reader.BCAFile()
    bca.open(filename)
    
    dss = bca_reader.DSS()
    alphabets = dss.get_alphabets()
    alphabet_sizes = [alphabet['size'] for alphabet in alphabets]
    
    # Compute the total size of the product alphabet
    product_alphabet_size = 1
    for size in alphabet_sizes:
        product_alphabet_size *= size
    print(f'{product_alphabet_size=}')
    
    # Initialize global counts counter (sparse representation)
    global_counts = Counter()
    
    for i in tqdm(range(len(bca)), total=len(bca), desc='Processing chains'):
        chain = bca[i]
        profiles = dss.get_profiles(chain)
        
        # Initialize accumulator with first profile as uint64
        accumulator = profiles[0].astype('uint64')
        
        # For each subsequent profile, multiply accumulator by next alphabet size and add profile
        for j in range(1, len(profiles)):
            accumulator = accumulator * alphabet_sizes[j] + profiles[j].astype('uint64')
        
        # Count occurrences of each alphabet symbol and add to global counter
        unique_values, counts = np.unique(accumulator, return_counts=True)
        for value, count in zip(unique_values, counts):
            global_counts[value] += count

    print('done accumulating')
    
    # Find the top k heavy hitters
    k = 1024
    top_k_items = global_counts.most_common(k)
    top_k_indices = np.array([item[0] for item in top_k_items])
    top_k_counts = np.array([item[1] for item in top_k_items])
    top_10_indices = top_k_indices[:10]  # Top 10 in descending order
    top_10_counts = top_k_counts[:10]
    
    print("Top 10 heavy hitters:")
    for rank, (idx, count) in enumerate(zip(top_10_indices, top_10_counts), 1):
        print(f"  Rank {rank}: symbol {idx}, count = {count}")
    
    # Calculate tail frequency (sum of frequencies not in top k)
    total_frequency = sum(global_counts.values())
    top_k_frequency = sum(top_k_counts)
    tail_frequency = total_frequency - top_k_frequency
    
    # Create rank-x semilogy plot for top k heavy hitters
    ranks = np.arange(1, k+1)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(ranks, top_k_counts, 'r-', linewidth=2)
    
    # Clean log scale with proper subdivisions
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.minorticks_on()
    ax.grid(True, which='minor', alpha=0.1)
    
    # Show minor tick labels for log scale
    ax.yaxis.set_minor_formatter(LogFormatter(minor_thresholds=(2, 0.4)))
    
    ax.set_xlabel('Rank')
    ax.set_ylabel('Frequency (log scale)')
    ax.set_title(f'Frequency vs Rank for Top {k} Heavy Hitters (Top-{k} freq: {top_k_frequency}, Tail freq: {tail_frequency})')
    plt.tight_layout()
    plt.savefig('heavy_hitters_plot.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'heavy_hitters_plot.png'")
    print(f"Total frequency: {total_frequency}")
    print(f"Top {k} frequency: {top_k_frequency}")
    print(f"Tail frequency: {tail_frequency}")
    
    bca.close()

if __name__ == "__main__":
    main()
