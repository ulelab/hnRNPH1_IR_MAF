#!/usr/bin/env python3
"""
Remove overlapping entries from a BED file.
For each gene (column 4), removes entries where start coordinates are within ±24 bp.
When duplicates are found, keeps the first occurrence.
"""

import pandas as pd
import argparse
import sys
import numpy as np

def remove_overlaps_by_gene(bed_df, window=24):
    """
    Remove entries where start coordinates are within ±window bp for the same gene.
    Groups by gene name and checks start coordinate proximity.
    """
    # Sort by gene, then by start position
    bed_df_sorted = bed_df.sort_values(['name', 'start', 'end']).reset_index(drop=True)
    
    # Track which rows to keep
    keep_mask = []
    
    # Process each gene group separately (much faster than checking all pairs)
    for gene_name, gene_group in bed_df_sorted.groupby('name'):
        gene_indices = gene_group.index.tolist()
        gene_starts = gene_group['start'].values
        
        # For this gene, mark which entries to keep
        gene_keep = [True] * len(gene_group)
        
        # Check each entry against previous entries in the same gene
        for i in range(1, len(gene_group)):
            current_start = gene_starts[i]
            
            # Check against all previous entries we're keeping
            for j in range(i):
                if not gene_keep[j]:
                    continue
                
                prev_start = gene_starts[j]
                
                # Check if start coordinates are within ±window bp
                if abs(current_start - prev_start) <= window:
                    gene_keep[i] = False
                    break  # Found a duplicate, no need to check further
        
        # Add to overall keep mask
        keep_mask.extend(gene_keep)
    
    # Convert to boolean array
    keep_mask = np.array(keep_mask)
    
    # Filter to keep only non-overlapping entries
    bed_df_no_overlaps = bed_df_sorted[keep_mask].copy()
    
    return bed_df_no_overlaps, sum(keep_mask), len(keep_mask) - sum(keep_mask)

def main():
    parser = argparse.ArgumentParser(
        description='Remove overlapping entries from a BED file'
    )
    parser.add_argument('--input', required=True,
                       help='Input BED file path')
    parser.add_argument('--output', required=True,
                       help='Output BED file path (without overlaps)')
    parser.add_argument('--window', type=int, default=24,
                       help='Window size in bp for detecting overlaps (default: 24)')
    parser.add_argument('--overlap_report', 
                       help='Optional: Output file to report overlapping entries')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Removing Overlapping Entries from BED File")
    print("=" * 60)
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print("=" * 60)
    print()
    
    # Read BED file
    print("Reading BED file...")
    try:
        bed_df = pd.read_csv(args.input, sep='\t', header=None,
                           names=['chr', 'start', 'end', 'name', 'score', 'strand'],
                           usecols=[0, 1, 2, 3, 4, 5])
    except Exception as e:
        print(f"ERROR reading BED file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"  Loaded {len(bed_df)} entries")
    print()
    
    # Remove overlaps
    print(f"Identifying and removing overlapping entries (window: ±{args.window} bp)...")
    bed_df_no_overlaps, n_kept, n_removed = remove_overlaps_by_gene(bed_df, window=args.window)
    
    print(f"  Entries kept: {n_kept}")
    print(f"  Entries removed: {n_removed}")
    print(f"  Overlap rate: {100 * n_removed / len(bed_df):.2f}%")
    print()
    
    # Write output
    print(f"Writing output to {args.output}...")
    bed_df_no_overlaps.to_csv(args.output, sep='\t', header=False, index=False)
    print(f"  Written {len(bed_df_no_overlaps)} entries")
    print()
    
    # Summary
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Input entries: {len(bed_df)}")
    print(f"Output entries: {len(bed_df_no_overlaps)}")
    print(f"Overlaps removed: {n_removed}")
    print(f"Output file: {args.output}")
    print("=" * 60)
    
    # If requested, report overlapping entries
    if args.overlap_report:
        print()
        print(f"Writing overlap report to {args.overlap_report}...")
        # Find which entries were removed (using same logic as removal)
        bed_df_sorted = bed_df.sort_values(['name', 'start', 'end']).reset_index(drop=True)
        
        removed_entries = []
        for gene_name, gene_group in bed_df_sorted.groupby('name'):
            gene_indices = gene_group.index.tolist()
            gene_starts = gene_group['start'].values
            gene_keep = [True] * len(gene_group)
            
            for i in range(1, len(gene_group)):
                current_start = gene_starts[i]
                current_row = gene_group.iloc[i]
                
                for j in range(i):
                    if not gene_keep[j]:
                        continue
                    
                    prev_start = gene_starts[j]
                    prev_row = gene_group.iloc[j]
                    
                    if abs(current_start - prev_start) <= args.window:
                        gene_keep[i] = False
                        removed_entries.append({
                            'chr': current_row['chr'],
                            'start': current_row['start'],
                            'end': current_row['end'],
                            'name': current_row['name'],
                            'score': current_row['score'],
                            'strand': current_row['strand'],
                            'overlaps_with_chr': prev_row['chr'],
                            'overlaps_with_start': prev_row['start'],
                            'overlaps_with_end': prev_row['end'],
                            'overlaps_with_name': prev_row['name'],
                            'distance_bp': abs(current_start - prev_start)
                        })
                        break
        
        if removed_entries:
            removed_df = pd.DataFrame(removed_entries)
            removed_df.to_csv(args.overlap_report, sep='\t', index=False)
            print(f"  Written {len(removed_entries)} overlapping entries to report")
        else:
            print("  No overlapping entries to report")
        print()

if __name__ == '__main__':
    main()

