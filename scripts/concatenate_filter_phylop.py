#!/usr/bin/env python3
"""
Concatenate all phyloP summary files and filter out exonic decoy IDs.
Preserves original number formatting.
"""

import pandas as pd
import sys
import glob
import os

def main():
    # File paths
    phylop_dir = "phyloP/decoy"
    overlap_file = "decoy_exon_overlaps.tsv"
    output_file = "phyloP_decoy_summary_all.tsv"
    filtered_output_file = "phyloP_decoy_summary_non_exonic.tsv"
    
    print("Loading overlap file to identify exonic decoy IDs...")
    # Read the overlap file to get decoy IDs that overlap with exons
    try:
        overlap_df = pd.read_csv(overlap_file, sep="\t")
        # Get unique decoy IDs that have exon overlaps
        overlapping_decoy_ids = set(overlap_df['decoy_name'].unique())
        print(f"  Found {len(overlapping_decoy_ids)} unique decoy IDs with exon overlaps")
    except FileNotFoundError:
        print(f"ERROR: Overlap file '{overlap_file}' not found.")
        print("  Please run the R Markdown file 'decoy_exon_overlaps.Rmd' first to generate it.")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read overlap file: {e}")
        sys.exit(1)
    
    print(f"\nFinding all phyloP summary files in {phylop_dir}...")
    # Find all phyloP summary files
    phylop_files = sorted(glob.glob(os.path.join(phylop_dir, "*_phylop_summary.tsv")))
    
    if not phylop_files:
        print(f"ERROR: No phyloP summary files found in {phylop_dir}")
        sys.exit(1)
    
    print(f"  Found {len(phylop_files)} phyloP summary files")
    
    print(f"\nConcatenating phyloP summary files...")
    # Read and concatenate all files
    all_dfs = []
    total_rows = 0
    
    for phylop_file in phylop_files:
        try:
            df = pd.read_csv(phylop_file, sep="\t", dtype=str)
            all_dfs.append(df)
            total_rows += len(df)
            print(f"  Loaded {os.path.basename(phylop_file)}: {len(df)} rows")
        except Exception as e:
            print(f"  WARNING: Failed to read {phylop_file}: {e}")
            continue
    
    if not all_dfs:
        print("ERROR: No phyloP files were successfully loaded")
        sys.exit(1)
    
    # Concatenate all dataframes
    combined_df = pd.concat(all_dfs, ignore_index=True)
    print(f"\n  Total rows after concatenation: {len(combined_df):,}")
    
    # Check column names
    print(f"  Columns: {', '.join(combined_df.columns[:5].tolist())}...")
    
    # Get the decoy ID column name (MAF_file in phyloP files)
    decoy_id_col = None
    for col in ['MAF_file', 'maf_file', 'decoyID', 'decoy_id', 'file']:
        if col in combined_df.columns:
            decoy_id_col = col
            break
    
    if decoy_id_col is None:
        # Assume first column is the decoy ID
        decoy_id_col = combined_df.columns[0]
        print(f"  Using first column '{decoy_id_col}' as decoy ID column")
    else:
        print(f"  Using column '{decoy_id_col}' as decoy ID column")
    
    # Extract decoy IDs from the column (remove .maf extension if present)
    if decoy_id_col in combined_df.columns:
        # Extract decoy ID from filename (remove .maf extension)
        combined_df['decoy_id_clean'] = combined_df[decoy_id_col].str.replace('.maf', '', regex=False).str.strip()
    else:
        print(f"ERROR: Decoy ID column '{decoy_id_col}' not found in data")
        sys.exit(1)
    
    print(f"\nSaving concatenated file: {output_file}")
    # Save concatenated file
    combined_df.to_csv(output_file, sep="\t", index=False, quoting=3)
    print(f"  Saved {len(combined_df):,} rows")
    
    print(f"\nFiltering out exonic decoy IDs...")
    # Filter to keep only non-exonic decoy IDs
    filtered_df = combined_df[~combined_df['decoy_id_clean'].isin(overlapping_decoy_ids)].copy()
    
    # Remove the temporary column
    filtered_df = filtered_df.drop(columns=['decoy_id_clean'])
    
    print(f"  Rows after filtering: {len(filtered_df):,}")
    print(f"  Rows removed: {len(combined_df) - len(filtered_df):,}")
    print(f"  Unique decoy IDs after filtering: {len(filtered_df[decoy_id_col].unique()):,}")
    
    print(f"\nSaving filtered file: {filtered_output_file}")
    # Save filtered file
    filtered_df.to_csv(filtered_output_file, sep="\t", index=False, quoting=3)
    
    print(f"  File saved successfully!")
    print(f"  Output file size: {len(filtered_df):,} rows x {len(filtered_df.columns)} columns")
    
    # Summary statistics
    print(f"\nSummary:")
    print(f"  Original decoy IDs: {len(combined_df[decoy_id_col].unique()):,}")
    print(f"  Exon-overlapping decoy IDs: {len(overlapping_decoy_ids):,}")
    print(f"  Non-exonic decoy IDs: {len(filtered_df[decoy_id_col].unique()):,}")
    print(f"  Percentage removed: {(len(overlapping_decoy_ids) / len(combined_df[decoy_id_col].unique()) * 100):.2f}%")

if __name__ == "__main__":
    main()
