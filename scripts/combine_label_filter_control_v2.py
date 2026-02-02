#!/usr/bin/env python3
"""
Combine control SpliceAI output TSV files, add origin suffixes (_C for canonical, _P for nonPRPF8),
and filter using decoy_exon_overlaps.tsv to remove exonic loci (n/2 from canonical, n/2 from nonPRPF8 per gene).
"""

import pandas as pd
import sys
import glob
import os
import random
from collections import Counter

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
    # File paths
    ctrl_output_dir = "spliceai_output/ctrl_output"
    canonical_bed = "Canonical_splice_sites.bed"
    nonprpf8_bed = "NonPRPF8SpliceSites_f250.bed"
    overlap_file = "decoy_exon_overlaps.tsv"
    output_combined = "control_spliceai_combined.tsv"
    output_filtered = "control_spliceai_filtered.tsv"
    
    # Set random seed for reproducibility
    random.seed(42)
    
    print("=" * 70)
    print("Combining and Filtering Control SpliceAI Dataset")
    print("=" * 70)
    print()
    
    # Step 1: Load source BED files to identify origin
    print("Step 1: Loading source BED files to identify origin...")
    canonical_ids = set()
    nonprpf8_ids = set()
    
    if os.path.exists(canonical_bed):
        try:
            try:
                canonical_df = pd.read_csv(canonical_bed, sep='\t', header=None,
                                          names=['chr', 'start', 'end', 'name', 'score', 'strand'],
                                          usecols=[0, 1, 2, 3, 4, 5])
            except:
                canonical_df = pd.read_csv(canonical_bed, sep=r'\s+', header=None,
                                          names=['chr', 'start', 'end', 'name', 'score', 'strand'],
                                          usecols=[0, 1, 2, 3, 4, 5])
            
            for _, row in canonical_df.iterrows():
                bed_name = str(row['name'])
                canonical_ids.add(bed_name)
                coord = int(row['start']) + (int(row['end']) - int(row['start'])) // 2
                decoy_id_center = f"{bed_name}_{coord}"
                canonical_ids.add(decoy_id_center)
                decoy_id_start = f"{bed_name}_{int(row['start'])}"
                canonical_ids.add(decoy_id_start)
            print(f"  Loaded {len(canonical_ids)} canonical decoy ID patterns")
        except Exception as e:
            print(f"  WARNING: Could not load canonical BED: {e}")
    else:
        print(f"  WARNING: Canonical BED file not found: {canonical_bed}")
    
    if os.path.exists(nonprpf8_bed):
        try:
            try:
                nonprpf8_df = pd.read_csv(nonprpf8_bed, sep='\t', header=None,
                                         names=['chr', 'start', 'end', 'name', 'score', 'strand'],
                                         usecols=[0, 1, 2, 3, 4, 5])
            except:
                nonprpf8_df = pd.read_csv(nonprpf8_bed, sep=r'\s+', header=None,
                                         names=['chr', 'start', 'end', 'name', 'score', 'strand'],
                                         usecols=[0, 1, 2, 3, 4, 5])
            
            for _, row in nonprpf8_df.iterrows():
                bed_name = str(row['name'])
                nonprpf8_ids.add(bed_name)
                coord = int(row['start']) + (int(row['end']) - int(row['start'])) // 2
                decoy_id_center = f"{bed_name}_{coord}"
                nonprpf8_ids.add(decoy_id_center)
                decoy_id_start = f"{bed_name}_{int(row['start'])}"
                nonprpf8_ids.add(decoy_id_start)
            print(f"  Loaded {len(nonprpf8_ids)} nonPRPF8 decoy ID patterns")
        except Exception as e:
            print(f"  WARNING: Could not load nonPRPF8 BED: {e}")
    else:
        print(f"  WARNING: NonPRPF8 BED file not found: {nonprpf8_bed}")
    
    print()
    
    # Step 2: Combine all control TSV files
    print("Step 2: Combining control SpliceAI TSV files...")
    tsv_files = sorted(glob.glob(os.path.join(ctrl_output_dir, "*_spliceai_scores.tsv")))
    
    if not tsv_files:
        print(f"ERROR: No TSV files found in {ctrl_output_dir}")
        sys.exit(1)
    
    print(f"  Found {len(tsv_files)} TSV files")
    
    all_dfs = []
    for tsv_file in tsv_files:
        try:
            df = pd.read_csv(tsv_file, sep='\t', dtype=str)
            all_dfs.append(df)
        except Exception as e:
            print(f"    WARNING: Failed to load {tsv_file}: {e}")
            continue
    
    if not all_dfs:
        print("ERROR: No TSV files were successfully loaded")
        sys.exit(1)
    
    combined_df = pd.concat(all_dfs, ignore_index=True)
    print(f"  Total rows after concatenation: {len(combined_df):,}")
    print()
    
    # Step 3: Add origin suffixes to decoyIDs
    print("Step 3: Adding origin suffixes to decoyIDs...")
    def add_suffix(decoy_id):
        """Add _C or _P suffix based on origin."""
        if decoy_id in canonical_ids:
            return f"{decoy_id}_C"
        elif decoy_id in nonprpf8_ids:
            return f"{decoy_id}_P"
        else:
            gene_name = parse_decoy_id(decoy_id)
            if gene_name:
                canonical_match = any(cid.startswith(f"{gene_name}_") for cid in canonical_ids)
                nonprpf8_match = any(pid.startswith(f"{gene_name}_") for pid in nonprpf8_ids)
                
                if canonical_match and nonprpf8_match:
                    return f"{decoy_id}_C"
                elif canonical_match:
                    return f"{decoy_id}_C"
                elif nonprpf8_match:
                    return f"{decoy_id}_P"
            return f"{decoy_id}_C"
    
    combined_df['decoyID'] = combined_df['decoyID'].apply(add_suffix)
    
    canonical_count = sum(1 for did in combined_df['decoyID'] if did.endswith('_C'))
    nonprpf8_count = sum(1 for did in combined_df['decoyID'] if did.endswith('_P'))
    print(f"  Canonical (_C): {canonical_count:,}")
    print(f"  NonPRPF8 (_P): {nonprpf8_count:,}")
    print()
    
    # Step 4: Format numeric columns
    print("Step 4: Formatting numeric columns...")
    numeric_cols = [col for col in combined_df.columns if col != 'decoyID']
    for col in numeric_cols:
        combined_df[col] = combined_df[col].apply(format_number)
    
    # Step 5: Save combined file
    print(f"Step 5: Saving combined file: {output_combined}")
    combined_df.to_csv(output_combined, sep='\t', index=False, quoting=3)
    print(f"  Saved {len(combined_df):,} rows")
    print()
    
    # Step 6: Load decoy-exon overlaps.tsv to identify genes with exonic loci
    print("Step 6: Loading decoy-exon overlaps...")
    
    # Count exonic loci per gene from the overlap file
    genes_with_exonic = Counter()  # {gene_name: count}
    
    if os.path.exists(overlap_file):
        try:
            overlap_df = pd.read_csv(overlap_file, sep='\t', dtype=str)
            # Get unique decoy names that overlap with exons
            exonic_decoy_ids = set(overlap_df['decoy_name'].unique())
            print(f"  Found {len(exonic_decoy_ids)} unique decoyIDs that overlap with exons")
            
            # Extract gene names from exonic decoyIDs
            for decoy_name in exonic_decoy_ids:
                gene_name = parse_decoy_id(decoy_name)
                if gene_name:
                    genes_with_exonic[gene_name] += 1
            
            print(f"  Found {len(genes_with_exonic)} genes with exonic loci")
            print(f"  Total exonic loci: {sum(genes_with_exonic.values())}")
        except Exception as e:
            print(f"ERROR: Failed to load overlap file: {e}")
            sys.exit(1)
    else:
        print(f"ERROR: Overlap file not found: {overlap_file}")
        sys.exit(1)
    
    print()
    
    # Step 7: Filter control dataset by removing exonic loci
    # For each gene with exonic loci, remove n/2 from canonical and n/2 from nonPRPF8
    print("Step 7: Filtering control dataset to remove exonic loci...")
    print("  Strategy: For each gene with N exonic loci, remove N/2 from canonical (_C) and N/2 from nonPRPF8 (_P)")
    print()
    
    # Extract gene names from control decoyIDs (remove suffix first)
    def get_gene_from_labeled_id(labeled_id):
        """Extract gene name from labeled decoyID (e.g., GENE_12345_C -> GENE)."""
        if labeled_id.endswith('_C'):
            base_id = labeled_id[:-2]
        elif labeled_id.endswith('_P'):
            base_id = labeled_id[:-2]
        else:
            base_id = labeled_id
        return parse_decoy_id(base_id)
    
    # Group control loci by gene and source type
    control_by_gene_source = {}  # {gene_name: {'_C': [decoy_ids], '_P': [decoy_ids]}}
    
    for labeled_id in combined_df['decoyID']:
        gene_name = get_gene_from_labeled_id(labeled_id)
        if not gene_name:
            continue
        
        if gene_name not in control_by_gene_source:
            control_by_gene_source[gene_name] = {'_C': [], '_P': []}
        
        if labeled_id.endswith('_C'):
            control_by_gene_source[gene_name]['_C'].append(labeled_id)
        elif labeled_id.endswith('_P'):
            control_by_gene_source[gene_name]['_P'].append(labeled_id)
    
    print(f"  Control dataset has loci from {len(control_by_gene_source)} genes")
    
    # For each gene with exonic loci, remove n/2 from canonical and n/2 from nonPRPF8
    decoy_ids_to_remove = set()
    genes_processed = 0
    genes_skipped = 0
    
    for gene_name, n_exonic in genes_with_exonic.items():
        if gene_name not in control_by_gene_source:
            genes_skipped += 1
            continue  # Gene not in control dataset
        
        genes_processed += 1
        canonical_loci = control_by_gene_source[gene_name]['_C']
        nonprpf8_loci = control_by_gene_source[gene_name]['_P']
        
        # Calculate how many to remove from each source (n/2 each)
        n_remove_canonical = n_exonic // 2
        n_remove_nonprpf8 = n_exonic // 2
        
        # If odd number, add the remainder to canonical
        if n_exonic % 2 == 1:
            n_remove_canonical += 1
        
        # Randomly sample which ones to remove
        if n_remove_canonical > 0 and len(canonical_loci) > 0:
            n_to_remove_c = min(n_remove_canonical, len(canonical_loci))
            to_remove_c = random.sample(canonical_loci, n_to_remove_c)
            decoy_ids_to_remove.update(to_remove_c)
        
        if n_remove_nonprpf8 > 0 and len(nonprpf8_loci) > 0:
            n_to_remove_p = min(n_remove_nonprpf8, len(nonprpf8_loci))
            to_remove_p = random.sample(nonprpf8_loci, n_to_remove_p)
            decoy_ids_to_remove.update(to_remove_p)
    
    print(f"  Genes processed: {genes_processed}")
    print(f"  Genes skipped (not in control): {genes_skipped}")
    print(f"  Total loci to remove: {len(decoy_ids_to_remove):,}")
    
    # Count by source
    removed_canonical = sum(1 for did in decoy_ids_to_remove if did.endswith('_C'))
    removed_nonprpf8 = sum(1 for did in decoy_ids_to_remove if did.endswith('_P'))
    print(f"    Canonical (_C): {removed_canonical:,}")
    print(f"    NonPRPF8 (_P): {removed_nonprpf8:,}")
    print()
    
    # Filter control dataset
    filtered_df = combined_df[~combined_df['decoyID'].isin(decoy_ids_to_remove)].copy()
    
    print(f"  Rows before filtering: {len(combined_df):,}")
    print(f"  Rows after filtering: {len(filtered_df):,}")
    print(f"  Rows removed: {len(combined_df) - len(filtered_df):,}")
    print()
    
    # Step 8: Save filtered file
    print(f"Step 8: Saving filtered file: {output_filtered}")
    filtered_df.to_csv(output_filtered, sep='\t', index=False, quoting=3)
    print(f"  Saved {len(filtered_df):,} rows")
    print()
    
    # Summary statistics
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Combined control dataset:")
    print(f"  Total rows: {len(combined_df):,}")
    print(f"  Canonical (_C): {canonical_count:,}")
    print(f"  NonPRPF8 (_P): {nonprpf8_count:,}")
    print()
    print(f"Exonic overlap filtering:")
    print(f"  Exonic decoyIDs found: {len(exonic_decoy_ids):,}")
    print(f"  Genes with exonic loci: {len(genes_with_exonic)}")
    print(f"  Genes processed: {genes_processed}")
    print(f"  Loci removed: {len(decoy_ids_to_remove):,}")
    print(f"    Canonical (_C): {removed_canonical:,}")
    print(f"    NonPRPF8 (_P): {removed_nonprpf8:,}")
    print()
    print(f"Filtered control dataset:")
    print(f"  Total rows: {len(filtered_df):,}")
    filtered_canonical = sum(1 for did in filtered_df['decoyID'] if did.endswith('_C'))
    filtered_nonprpf8 = sum(1 for did in filtered_df['decoyID'] if did.endswith('_P'))
    print(f"  Canonical (_C): {filtered_canonical:,}")
    print(f"  NonPRPF8 (_P): {filtered_nonprpf8:,}")
    print()
    print(f"Output files:")
    print(f"  Combined: {output_combined}")
    print(f"  Filtered: {output_filtered}")
    print("=" * 70)

if __name__ == "__main__":
    main()
