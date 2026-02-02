#!/usr/bin/env python3
"""
Create a control dataset that matches gene representation in the experimental dataset.
1. Count gene occurrences in SpliceAI_MAF.tsv
2. For each gene, sample half (rounded up) occurrences from Canonical_splice_sites.bed
3. For each gene, sample half (rounded up) occurrences from NonPRPF8SpliceSites_f250.bed
4. Combine, remove duplicates, expand coordinates, and output
"""

import pandas as pd
import numpy as np
import argparse
import sys
import math
from collections import Counter

def parse_decoy_id(decoy_id):
    """
    Parse decoy ID like 'A1BG-AS1_58353595' into (gene_name, coordinate).
    The gene name is everything before the last underscore.
    """
    parts = decoy_id.rsplit('_', 1)
    if len(parts) != 2:
        return None
    gene_name = parts[0]
    return gene_name

def expand_coordinates(df, expand_bp=24):
    """
    Expand BED coordinates by expand_bp on each side.
    Ensures start doesn't go below 0.
    """
    df_expanded = df.copy()
    df_expanded['start'] = (df_expanded['start'] - expand_bp).clip(lower=0)
    df_expanded['end'] = df_expanded['end'] + expand_bp
    return df_expanded

def main():
    parser = argparse.ArgumentParser(
        description='Create control dataset matching gene representation in experimental data'
    )
    parser.add_argument('--experimental_maf', required=True,
                       help='Path to SpliceAI_MAF.tsv (experimental dataset)')
    parser.add_argument('--canonical', required=True,
                       help='Path to Canonical_splice_sites.bed')
    parser.add_argument('--nonprpf8', required=True,
                       help='Path to NonPRPF8SpliceSites_f250.bed')
    parser.add_argument('--output', required=True,
                       help='Output BED file path')
    parser.add_argument('--expand', type=int, default=24,
                       help='Base pairs to expand on each side (default: 24)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--min_score', type=float, default=0.0002,
                       help='Minimum score for non-PRPF8 sites (default: 0.0002)')
    
    args = parser.parse_args()
    
    # Set random seed
    np.random.seed(args.seed)
    
    print("=" * 60)
    print("Creating Matched Control Dataset")
    print("=" * 60)
    print(f"Experimental MAF: {args.experimental_maf}")
    print(f"Canonical file: {args.canonical}")
    print(f"Non-PRPF8 file: {args.nonprpf8}")
    print(f"Coordinate expansion: ±{args.expand} bp")
    print("=" * 60)
    print()
    
    # Step 1: Count gene occurrences in experimental dataset
    print("Step 1: Counting gene occurrences in experimental dataset...")
    try:
        experimental_df = pd.read_csv(args.experimental_maf, sep='\t', header=0)
    except Exception as e:
        print(f"ERROR reading experimental MAF: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"  Loaded {len(experimental_df)} decoy IDs from experimental dataset")
    
    # Parse gene names from decoy IDs
    gene_counts = Counter()
    for decoy_id in experimental_df['decoyID']:
        gene_name = parse_decoy_id(decoy_id)
        if gene_name:
            gene_counts[gene_name] += 1
    
    print(f"  Found {len(gene_counts)} unique genes")
    print(f"  Total gene occurrences: {sum(gene_counts.values())}")
    print()
    
    # Step 2: Load canonical and non-PRPF8 BED files
    print("Step 2: Loading source BED files...")
    
    print(f"Reading {args.canonical}...")
    try:
        canonical = pd.read_csv(args.canonical, sep='\t', header=None,
                               names=['chr', 'start', 'end', 'name', 'score', 'strand'],
                               usecols=[0, 1, 2, 3, 4, 5])
    except Exception as e:
        print(f"ERROR reading canonical file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"  Loaded {len(canonical)} rows")
    
    print(f"Reading {args.nonprpf8}...")
    try:
        nonprpf8 = pd.read_csv(args.nonprpf8, sep='\t', header=None,
                               names=['chr', 'start', 'end', 'name', 'score', 'strand'],
                               usecols=[0, 1, 2, 3, 4, 5])
    except Exception as e:
        print(f"ERROR reading non-PRPF8 file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"  Loaded {len(nonprpf8)} rows")
    
    # Filter non-PRPF8 sites: score > min_score
    print(f"Filtering non-PRPF8 sites: score > {args.min_score}...")
    nonprpf8_filtered = nonprpf8[nonprpf8['score'] > args.min_score].copy()
    print(f"  After filtering: {len(nonprpf8_filtered)} rows (removed {len(nonprpf8) - len(nonprpf8_filtered)})")
    nonprpf8 = nonprpf8_filtered
    print()
    
    # Step 3: Sample from each source based on gene counts
    print("Step 3: Sampling from source files based on gene representation...")
    
    canonical_samples = []
    nonprpf8_samples = []
    
    genes_with_canonical = 0
    genes_with_nonprpf8 = 0
    genes_with_both = 0
    genes_with_neither = 0
    
    total_canonical_requested = 0
    total_nonprpf8_requested = 0
    total_canonical_sampled = 0
    total_nonprpf8_sampled = 0
    
    for gene_name, count in gene_counts.items():
        # Calculate how many to sample: half rounded up
        n_sample = math.ceil(count / 2)
        
        # Get entries for this gene from canonical
        canonical_gene = canonical[canonical['name'] == gene_name]
        # Get entries for this gene from non-PRPF8
        nonprpf8_gene = nonprpf8[nonprpf8['name'] == gene_name]
        
        has_canonical = len(canonical_gene) > 0
        has_nonprpf8 = len(nonprpf8_gene) > 0
        
        if has_canonical and has_nonprpf8:
            genes_with_both += 1
        elif has_canonical:
            genes_with_canonical += 1
        elif has_nonprpf8:
            genes_with_nonprpf8 += 1
        else:
            genes_with_neither += 1
        
        # Sample from canonical
        if has_canonical:
            total_canonical_requested += n_sample
            n_canonical_sample = min(n_sample, len(canonical_gene))
            if n_canonical_sample > 0:
                canonical_sample = canonical_gene.sample(n=n_canonical_sample, random_state=args.seed)
                canonical_samples.append(canonical_sample)
                total_canonical_sampled += n_canonical_sample
        
        # Sample from non-PRPF8
        if has_nonprpf8:
            total_nonprpf8_requested += n_sample
            n_nonprpf8_sample = min(n_sample, len(nonprpf8_gene))
            if n_nonprpf8_sample > 0:
                nonprpf8_sample = nonprpf8_gene.sample(n=n_nonprpf8_sample, random_state=args.seed + 1)
                nonprpf8_samples.append(nonprpf8_sample)
                total_nonprpf8_sampled += n_nonprpf8_sample
    
    print(f"  Genes with canonical sites: {genes_with_canonical}")
    print(f"  Genes with non-PRPF8 sites: {genes_with_nonprpf8}")
    print(f"  Genes with both: {genes_with_both}")
    print(f"  Genes with neither: {genes_with_neither}")
    print()
    print(f"  Canonical: requested {total_canonical_requested}, sampled {total_canonical_sampled}")
    print(f"  Non-PRPF8: requested {total_nonprpf8_requested}, sampled {total_nonprpf8_sampled}")
    print()
    
    # Combine samples
    print("Step 4: Combining samples...")
    if canonical_samples:
        canonical_combined = pd.concat(canonical_samples, ignore_index=True)
    else:
        canonical_combined = pd.DataFrame(columns=['chr', 'start', 'end', 'name', 'score', 'strand'])
    
    if nonprpf8_samples:
        nonprpf8_combined = pd.concat(nonprpf8_samples, ignore_index=True)
    else:
        nonprpf8_combined = pd.DataFrame(columns=['chr', 'start', 'end', 'name', 'score', 'strand'])
    
    combined = pd.concat([canonical_combined, nonprpf8_combined], ignore_index=True)
    print(f"  Combined total: {len(combined)} rows")
    print()
    
    # Check for duplicates (by chr, start, end, strand)
    print("Step 5: Checking for duplicates...")
    duplicate_key = ['chr', 'start', 'end', 'strand']
    duplicates = combined.duplicated(subset=duplicate_key, keep=False)
    n_duplicates = duplicates.sum()
    
    if n_duplicates > 0:
        print(f"  Found {n_duplicates} duplicate rows (same chr, start, end, strand)")
        print("  Removing duplicates (keeping first occurrence)...")
        combined = combined.drop_duplicates(subset=duplicate_key, keep='first')
        print(f"  After deduplication: {len(combined)} rows")
    else:
        print("  No duplicates found")
    print()
    
    # Expand coordinates
    print(f"Step 6: Expanding coordinates by ±{args.expand} bp...")
    combined_expanded = expand_coordinates(combined, expand_bp=args.expand)
    
    # Show some statistics
    print()
    print("Coordinate expansion summary:")
    original_widths = combined['end'] - combined['start']
    expanded_widths = combined_expanded['end'] - combined_expanded['start']
    print(f"  Original width range: {original_widths.min()} - {original_widths.max()} bp")
    print(f"  Expanded width range: {expanded_widths.min()} - {expanded_widths.max()} bp")
    print(f"  Mean expansion: {(expanded_widths - original_widths).mean():.1f} bp")
    print()
    
    # Sort by chromosome and start position
    print("Step 7: Sorting by chromosome and start position...")
    combined_expanded = combined_expanded.sort_values(['chr', 'start', 'end'])
    print()
    
    # Write output
    print(f"Step 8: Writing output to {args.output}...")
    combined_expanded.to_csv(args.output, sep='\t', header=False, index=False)
    print(f"  Written {len(combined_expanded)} rows")
    print()
    
    # Summary
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Experimental dataset genes: {len(gene_counts)}")
    print(f"Experimental dataset total occurrences: {sum(gene_counts.values())}")
    print()
    print(f"Canonical samples: {len(canonical_combined)}")
    print(f"Non-PRPF8 samples: {len(nonprpf8_combined)}")
    print(f"Combined (before dedup): {len(canonical_combined) + len(nonprpf8_combined)}")
    print(f"Final (after dedup): {len(combined_expanded)}")
    print(f"Duplicates removed: {len(canonical_combined) + len(nonprpf8_combined) - len(combined_expanded)}")
    print()
    print(f"Gene coverage:")
    print(f"  Genes with canonical sites only: {genes_with_canonical}")
    print(f"  Genes with non-PRPF8 sites only: {genes_with_nonprpf8}")
    print(f"  Genes with both: {genes_with_both}")
    print(f"  Genes with neither: {genes_with_neither}")
    print()
    print(f"Output file: {args.output}")
    print("=" * 60)

if __name__ == '__main__':
    main()

