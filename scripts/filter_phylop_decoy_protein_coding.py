#!/usr/bin/env python3
"""
Filter phyloP_decoy_summary_all.tsv to match decoyIDs in SpliceAI_MAF_protein_coding.tsv
"""

import pandas as pd

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
    phylop_file = "phyloP/decoy/phyloP_decoy_summary_all.tsv"
    spliceai_file = "SpliceAI_MAF_protein_coding.tsv"
    output_file = "phyloP/decoy/phyloP_decoy_protein_coding.tsv"
    
    print("Loading SpliceAI protein coding file to get decoy IDs...")
    spliceai_df = pd.read_csv(spliceai_file, sep="\t", dtype=str, usecols=['decoyID'])
    spliceai_decoy_ids = set(spliceai_df['decoyID'].unique())
    print(f"  Found {len(spliceai_decoy_ids)} unique decoy IDs in SpliceAI file")
    
    print(f"\nLoading phyloP decoy summary file...")
    phylop_df = pd.read_csv(phylop_file, sep="\t", dtype=str)
    print(f"  Total rows: {len(phylop_df):,}")
    print(f"  Columns: {', '.join(phylop_df.columns.tolist())}")
    
    # Get the decoy ID column name (MAF_file in phyloP files)
    decoy_id_col = 'MAF_file' if 'MAF_file' in phylop_df.columns else phylop_df.columns[0]
    print(f"  Using column '{decoy_id_col}' as decoy ID column")
    
    # Filter to keep only decoy IDs that are in the SpliceAI file
    print(f"\nFiltering to match SpliceAI decoy IDs...")
    filtered_df = phylop_df[phylop_df[decoy_id_col].isin(spliceai_decoy_ids)].copy()
    
    print(f"  Rows after filtering: {len(filtered_df):,}")
    print(f"  Rows removed: {len(phylop_df) - len(filtered_df):,}")
    print(f"  Unique decoy IDs after filtering: {len(filtered_df[decoy_id_col].unique()):,}")
    
    # Format numeric columns
    numeric_cols = ['mean_phylop', 'max_phylop', 'fraction_constrained', 
                    'fraction_constrained_p01', 'sum_positive_phylop']
    
    for col in numeric_cols:
        if col in filtered_df.columns:
            filtered_df[col] = filtered_df[col].apply(format_number)
    
    print(f"\nSaving filtered file: {output_file}")
    filtered_df.to_csv(output_file, sep="\t", index=False, quoting=3)
    
    print(f"  File saved successfully!")
    print(f"  Output file size: {len(filtered_df):,} rows x {len(filtered_df.columns)} columns")

if __name__ == "__main__":
    main()
