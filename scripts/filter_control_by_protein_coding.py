#!/usr/bin/env python3
"""
Filter control_spliceai_combined.tsv to keep only decoyIDs that match gene representation in
SpliceAI_MAF_protein_coding.tsv (matching by gene name, not exact decoyID).
"""

import pandas as pd
import sys

def parse_decoy_id(decoy_id):
    """Parse decoy ID to extract gene name (everything before last underscore)."""
    parts = decoy_id.rsplit('_', 1)
    if len(parts) != 2:
        return None
    return parts[0]

def format_number(value):
    """Format number to preserve decimal format (no scientific notation)."""
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
    input_control = "control_spliceai_combined.tsv"
    filter_file = "SpliceAI_MAF_protein_coding.tsv"
    output_filtered = "control_spliceai_filtered_protein_coding.tsv"
    
    print("=" * 70)
    print("Filtering Control Dataset by Protein Coding Gene Representation")
    print("=" * 70)
    print()
    
    # Step 1: Load filter file (protein coding)
    print("Step 1: Loading filter file...")
    try:
        filter_df = pd.read_csv(filter_file, sep='\t', dtype=str)
        filter_decoy_ids = set(filter_df['decoyID'].unique())
        print(f"  Loaded {len(filter_df):,} rows from {filter_file}")
        print(f"  Found {len(filter_decoy_ids):,} unique decoy IDs in filter file")
    except Exception as e:
        print(f"ERROR: Failed to load filter file: {e}")
        sys.exit(1)
    print()
    
    # Step 2: Load control dataset
    print("Step 2: Loading control dataset...")
    try:
        control_df = pd.read_csv(input_control, sep='\t', dtype=str)
        print(f"  Loaded {len(control_df):,} rows from {input_control}")
    except Exception as e:
        print(f"ERROR: Failed to load control file: {e}")
        sys.exit(1)
    print()
    
    # Step 3: Filter control dataset by gene name
    print("Step 3: Filtering control dataset...")
    print("  Matching by gene name (extracted from decoyID)...")
    
    # Extract gene names from filter file
    filter_genes = set()
    for decoy_id in filter_decoy_ids:
        gene_name = parse_decoy_id(decoy_id)
        if gene_name:
            filter_genes.add(gene_name)
    
    print(f"  Filter file has {len(filter_genes)} unique genes")
    
    # Extract gene names from control dataset
    def get_gene_from_labeled_id(labeled_id):
        """Extract gene name from labeled decoyID (e.g., GENE_12345_C -> GENE)."""
        if labeled_id.endswith('_C'):
            base_id = labeled_id[:-2]
        elif labeled_id.endswith('_P'):
            base_id = labeled_id[:-2]
        else:
            base_id = labeled_id
        return parse_decoy_id(base_id)
    
    control_genes = {}
    for labeled_id in control_df['decoyID']:
        gene_name = get_gene_from_labeled_id(labeled_id)
        if gene_name:
            if gene_name not in control_genes:
                control_genes[gene_name] = []
            control_genes[gene_name].append(labeled_id)
    
    print(f"  Control dataset has {len(control_genes)} unique genes")
    
    # Find which genes from control are in the filter file
    matching_genes = set(control_genes.keys()) & filter_genes
    print(f"  Found {len(matching_genes)} matching genes")
    
    # Get all labeled IDs that match
    matching_labeled_ids = set()
    for gene_name in matching_genes:
        matching_labeled_ids.update(control_genes[gene_name])
    
    print(f"  Found {len(matching_labeled_ids)} labeled decoyIDs (with _C/_P) to keep")
    
    # Filter the dataframe
    filtered_df = control_df[control_df['decoyID'].isin(matching_labeled_ids)].copy()
    
    print(f"  Filtered dataset: {len(filtered_df):,} rows")
    print()
    
    # Count by suffix
    canonical_count = sum(filtered_df['decoyID'].str.endswith('_C'))
    nonprpf8_count = sum(filtered_df['decoyID'].str.endswith('_P'))
    print(f"  Canonical (_C): {canonical_count:,} rows")
    print(f"  NonPRPF8 (_P): {nonprpf8_count:,} rows")
    print()
    
    # Step 4: Format numeric columns
    print("Step 4: Formatting numeric columns...")
    numeric_cols = [col for col in filtered_df.columns if col != 'decoyID']
    for col in numeric_cols:
        filtered_df[col] = filtered_df[col].apply(format_number)
    print("  Formatted numeric columns to fixed-point decimal format")
    print()
    
    # Step 5: Save filtered file
    print("Step 5: Saving filtered file...")
    try:
        filtered_df.to_csv(output_filtered, sep='\t', index=False, quoting=3)
        print(f"  Saved {len(filtered_df):,} rows to {output_filtered}")
    except Exception as e:
        print(f"ERROR: Failed to save filtered file: {e}")
        sys.exit(1)
    print()
    
    print("=" * 70)
    print("SUCCESS: Filtering completed")
    print("=" * 70)
    print(f"Output file: {output_filtered}")
    print(f"  Total rows: {len(filtered_df):,}")
    print(f"  Canonical (_C): {canonical_count:,}")
    print(f"  NonPRPF8 (_P): {nonprpf8_count:,}")
    print()

if __name__ == "__main__":
    main()
