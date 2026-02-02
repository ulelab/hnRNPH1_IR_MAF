#!/usr/bin/env python3
"""
Split SpliceAI_MAF_non_exonic.tsv into two files based on unstranded exon overlaps:
1. decoyIDs that have unstranded overlaps with exons
2. decoyIDs that do not have unstranded overlaps with exons
"""

import pandas as pd
import sys

def format_number(value):
    """Format numeric values to fixed-point decimal (6 decimal places) to avoid scientific notation."""
    if isinstance(value, str):
        if value == "-1.000000" or value == "-1.0" or value == "-1":
            return "-1.000000"
        try:
            float_val = float(value)
            if float_val == -1.0:
                return "-1.000000"
            return f"{float_val:.6f}"
        except (ValueError, TypeError):
            return value
    if pd.isna(value):
        return "-1.000000"
    try:
        float_val = float(value)
        if float_val == -1.0:
            return "-1.000000"
        return f"{float_val:.6f}"
    except (ValueError, TypeError):
        return str(value)

def main():
    overlap_file = "decoy_exon_overlaps_unstranded.tsv"
    input_file = "SpliceAI_MAF_non_exonic.tsv"
    output_overlaps = "SpliceAI_MAF_non_exonic_unstranded_overlaps.tsv"
    output_no_overlaps = "SpliceAI_MAF_non_exonic_no_unstranded_overlaps.tsv"
    
    print("=" * 70)
    print("Splitting SpliceAI_MAF_non_exonic.tsv by unstranded exon overlaps")
    print("=" * 70)
    print()
    
    # Step 1: Load unstranded overlaps file
    print("Step 1: Loading unstranded overlaps file...")
    try:
        overlap_df = pd.read_csv(overlap_file, sep='\t', dtype=str)
        overlapping_decoy_ids = set(overlap_df['decoy_name'].unique())
        print(f"  Loaded {len(overlap_df):,} overlap records")
        print(f"  Found {len(overlapping_decoy_ids):,} unique decoy IDs with unstranded overlaps")
    except Exception as e:
        print(f"ERROR: Failed to load overlap file: {e}")
        sys.exit(1)
    print()
    
    # Step 2: Load SpliceAI_MAF_non_exonic.tsv
    print("Step 2: Loading SpliceAI_MAF_non_exonic.tsv...")
    try:
        # Read with dtype=str to preserve original formatting
        spliceai_df = pd.read_csv(input_file, sep='\t', dtype=str)
        print(f"  Loaded {len(spliceai_df):,} rows")
        print(f"  Columns: {', '.join(spliceai_df.columns[:5])}...")
    except Exception as e:
        print(f"ERROR: Failed to load input file: {e}")
        sys.exit(1)
    print()
    
    # Step 3: Split into overlapping and non-overlapping
    print("Step 3: Splitting data by unstranded overlaps...")
    overlapping_mask = spliceai_df['decoyID'].isin(overlapping_decoy_ids)
    overlapping_df = spliceai_df[overlapping_mask].copy()
    no_overlaps_df = spliceai_df[~overlapping_mask].copy()
    
    print(f"  Decoy IDs with unstranded overlaps: {len(overlapping_df):,}")
    print(f"  Decoy IDs without unstranded overlaps: {len(no_overlaps_df):,}")
    print(f"  Total: {len(overlapping_df) + len(no_overlaps_df):,}")
    print()
    
    # Step 4: Format numeric columns
    print("Step 4: Formatting numeric columns...")
    numeric_cols = [col for col in overlapping_df.columns if col != 'decoyID']
    
    for col in numeric_cols:
        overlapping_df[col] = overlapping_df[col].apply(format_number)
        no_overlaps_df[col] = no_overlaps_df[col].apply(format_number)
    
    print("  Formatted numeric columns to fixed-point decimal format")
    print()
    
    # Step 5: Save output files
    print("Step 5: Saving output files...")
    try:
        overlapping_df.to_csv(output_overlaps, sep='\t', index=False, quoting=3)
        print(f"  Saved {len(overlapping_df):,} rows to {output_overlaps}")
    except Exception as e:
        print(f"ERROR: Failed to save overlapping file: {e}")
        sys.exit(1)
    
    try:
        no_overlaps_df.to_csv(output_no_overlaps, sep='\t', index=False, quoting=3)
        print(f"  Saved {len(no_overlaps_df):,} rows to {output_no_overlaps}")
    except Exception as e:
        print(f"ERROR: Failed to save non-overlapping file: {e}")
        sys.exit(1)
    print()
    
    print("=" * 70)
    print("SUCCESS: Split completed")
    print("=" * 70)
    print(f"Output files:")
    print(f"  - {output_overlaps}: {len(overlapping_df):,} rows (with unstranded overlaps)")
    print(f"  - {output_no_overlaps}: {len(no_overlaps_df):,} rows (without unstranded overlaps)")
    print()

if __name__ == "__main__":
    main()
