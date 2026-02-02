#!/usr/bin/env python3
"""
Concatenate all phyloP control summary files, add _C and _P suffixes,
and filter to match IDs in control_spliceai_filtered_protein_coding.tsv
"""

import pandas as pd
import glob
import os

def format_number(value):
    """Format number to fixed-point decimal, handling -1.0 and scientific notation"""
    try:
        num = float(value)
        if num == -1.0:
            return "-1.000000"
        # Format to 6 decimal places
        return f"{num:.6f}"
    except (ValueError, TypeError):
        return str(value)

def main():
    phylop_dir = "phyloP/ctrl"
    control_file = "control_spliceai_filtered_protein_coding.tsv"
    output_file = "phyloP_ctrl_summary_filtered_protein_coding.tsv"
    
    print("Loading control SpliceAI file to get decoy IDs...")
    control_df = pd.read_csv(control_file, sep="\t", dtype=str)
    control_decoy_ids = set(control_df['decoyID'].unique())
    print(f"  Found {len(control_decoy_ids)} unique decoy IDs in control file")
    
    # Create a mapping from base decoyID to suffix (_C or _P)
    base_to_suffix = {}
    for decoy_id in control_decoy_ids:
        if decoy_id.endswith('_C'):
            base = decoy_id[:-2]
            base_to_suffix[base] = '_C'
        elif decoy_id.endswith('_P'):
            base = decoy_id[:-2]
            base_to_suffix[base] = '_P'
    
    print(f"  Found {len(base_to_suffix)} base decoy IDs with suffixes")
    
    print(f"\nFinding all phyloP control summary files in {phylop_dir}...")
    phylop_files = sorted(glob.glob(os.path.join(phylop_dir, "*_phylop_summary.tsv")))
    
    if not phylop_files:
        print(f"ERROR: No phyloP summary files found in {phylop_dir}")
        return
    
    print(f"  Found {len(phylop_files)} phyloP summary files")
    
    print(f"\nConcatenating phyloP control summary files...")
    all_dfs = []
    
    for phylop_file in phylop_files:
        try:
            df = pd.read_csv(phylop_file, sep="\t", dtype=str)
            all_dfs.append(df)
            print(f"  Loaded {os.path.basename(phylop_file)}: {len(df)} rows")
        except Exception as e:
            print(f"  WARNING: Failed to read {phylop_file}: {e}")
            continue
    
    if not all_dfs:
        print("ERROR: No phyloP files were successfully loaded")
        return
    
    # Concatenate all dataframes
    combined_df = pd.concat(all_dfs, ignore_index=True)
    print(f"\n  Total rows after concatenation: {len(combined_df):,}")
    
    # Get the decoy ID column name (MAF_file in phyloP files)
    decoy_id_col = 'MAF_file' if 'MAF_file' in combined_df.columns else combined_df.columns[0]
    print(f"  Using column '{decoy_id_col}' as decoy ID column")
    
    # Add suffixes to decoy IDs
    print(f"\nAdding _C and _P suffixes to decoy IDs...")
    def add_suffix(decoy_id):
        base = str(decoy_id).strip()
        if base in base_to_suffix:
            return base + base_to_suffix[base]
        return base
    
    combined_df['decoyID'] = combined_df[decoy_id_col].apply(add_suffix)
    
    # Filter to keep only decoy IDs that are in the control file
    print(f"\nFiltering to match control decoy IDs...")
    filtered_df = combined_df[combined_df['decoyID'].isin(control_decoy_ids)].copy()
    
    # Drop the original MAF_file column if it exists and is different from decoyID
    if decoy_id_col != 'decoyID' and decoy_id_col in filtered_df.columns:
        filtered_df = filtered_df.drop(columns=[decoy_id_col])
    
    print(f"  Rows after filtering: {len(filtered_df):,}")
    print(f"  Rows removed: {len(combined_df) - len(filtered_df):,}")
    print(f"  Unique decoy IDs after filtering: {len(filtered_df['decoyID'].unique()):,}")
    
    # Format numeric columns
    numeric_cols = ['mean_phylop', 'max_phylop', 'fraction_constrained', 
                    'fraction_constrained_p01', 'sum_positive_phylop']
    
    for col in numeric_cols:
        if col in filtered_df.columns:
            filtered_df[col] = filtered_df[col].apply(format_number)
    
    # Reorder columns to put decoyID first
    cols = ['decoyID'] + [c for c in filtered_df.columns if c != 'decoyID']
    filtered_df = filtered_df[cols]
    
    print(f"\nSaving filtered file: {output_file}")
    filtered_df.to_csv(output_file, sep="\t", index=False, quoting=3)
    
    print(f"  File saved successfully!")
    print(f"  Output file size: {len(filtered_df):,} rows x {len(filtered_df.columns)} columns")
    
    # Summary by suffix
    canonical_count = sum(filtered_df['decoyID'].str.endswith('_C'))
    nonprpf8_count = sum(filtered_df['decoyID'].str.endswith('_P'))
    print(f"\nSummary by origin:")
    print(f"  Canonical (_C): {canonical_count}")
    print(f"  NonPRPF8 (_P): {nonprpf8_count}")

if __name__ == "__main__":
    main()
