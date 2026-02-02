#!/usr/bin/env python3
"""
Filter SpliceAI_MAF.tsv to keep rows where value is either -1 or > 0.05 
in 8 or more of specified fish species, and replace species IDs with common names.
"""

import pandas as pd
import sys
import argparse

# Target species (common names -> species IDs)
TARGET_SPECIES = {
    'Lamprey': 'petMar2',
    'Spotted gar': 'lepOcu1',
    'Zebrafish': 'danRer10',
    'Mexican tetra': 'astMex1',
    'Atlantic cod': 'gadMor1',
    'Stickleback': 'gasAcu1',
    'Medaka': 'oryLat2',
    'Southern platyfish': 'xipMac1',
    'Tetraodon': 'tetNig2',
    'Fugu': 'fr3',
    'Yellowbelly pufferfish': 'takFla1',
    'Nile tilapia': 'oreNil2',
    "Princess of Burundi": 'neoBri1',
    "Burton's mouthbreeder": 'hapBur1',
    'Zebra mbuna': 'mayZeb1',
    'Pundamilia nyererei': 'punNye1',
    'Coelacanth': 'latCha1'
}

def load_species_mapping(mapping_file):
    """Load species ID to common name mapping."""
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    # Create mapping from species_id to matched_common_name
    id_to_name = dict(zip(mapping_df['species_id'], mapping_df['matched_common_name']))
    return id_to_name

def filter_row(row, target_species_ids, min_species=8, min_conserved=15, max_missing=3, threshold=0.1):
    """
    Check if row meets conservation criteria:
    1. Has value > threshold in at least min_conserved species (conserved sites)
    2. Has value = -1 in fewer than max_missing+1 species (not too many missing)
    
    Returns True if both conditions are met.
    """
    conserved_count = 0
    missing_count = 0
    
    for species_id in target_species_ids:
        if species_id in row.index:
            value = row[species_id]
            
            # Count conserved sites (value > threshold)
            if isinstance(value, (int, float)) and value > threshold:
                conserved_count += 1
            
            # Count missing sites (value = -1)
            if value == -1.0:
                missing_count += 1
    
    # Must have at least min_conserved conserved sites AND fewer than max_missing+1 missing sites
    return conserved_count >= min_conserved and missing_count <= max_missing

def main():
    parser = argparse.ArgumentParser(
        description='Filter SpliceAI_MAF.tsv by conservation criteria and replace species IDs with common names'
    )
    parser.add_argument('--input', required=True,
                       help='Input SpliceAI_MAF.tsv file')
    parser.add_argument('--mapping', required=True,
                       help='tree_to_clade_mapping.tsv file')
    parser.add_argument('--output', required=True,
                       help='Output filtered TSV file')
    parser.add_argument('--min-conserved', type=int, default=15,
                       help='Minimum number of species with value > threshold (default: 15)')
    parser.add_argument('--max-missing', type=int, default=3,
                       help='Maximum number of species with value = -1 (default: 3, i.e., < 4)')
    parser.add_argument('--threshold', type=float, default=0.1,
                       help='SpliceAI threshold value (default: 0.1)')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("Filtering Conservation Data")
    print("=" * 70)
    print(f"Input file: {args.input}")
    print(f"Mapping file: {args.mapping}")
    print(f"Output file: {args.output}")
    print(f"Minimum conserved species (value > {args.threshold}): {args.min_conserved}")
    print(f"Maximum missing species (value = -1): {args.max_missing}")
    print(f"SpliceAI threshold: {args.threshold}")
    print("=" * 70)
    print()
    
    # Load species mapping
    print("Loading species mapping...")
    id_to_name = load_species_mapping(args.mapping)
    print(f"  Loaded {len(id_to_name)} species mappings")
    print()
    
    # Get target species IDs
    target_species_ids = list(TARGET_SPECIES.values())
    target_species_names = list(TARGET_SPECIES.keys())
    
    print(f"Target species ({len(target_species_ids)}):")
    for name, sid in TARGET_SPECIES.items():
        print(f"  {name}: {sid}")
    print()
    
    # Load SpliceAI MAF data
    print("Loading SpliceAI MAF data...")
    df = pd.read_csv(args.input, sep='\t')
    print(f"  Loaded {len(df)} rows, {len(df.columns)} columns")
    print(f"  Columns: {list(df.columns[:5])}... (showing first 5)")
    print()
    
    # Check which target species are present in the data
    missing_species = []
    present_species = []
    for sid in target_species_ids:
        if sid in df.columns:
            present_species.append(sid)
        else:
            missing_species.append(sid)
    
    if missing_species:
        print(f"WARNING: {len(missing_species)} target species not found in data:")
        for sid in missing_species:
            name = [k for k, v in TARGET_SPECIES.items() if v == sid][0]
            print(f"  {name} ({sid})")
        print()
    
    print(f"Found {len(present_species)} target species in data")
    print()
    
    # Filter rows
    print(f"Filtering rows:")
    print(f"  - Value > {args.threshold} in >= {args.min_conserved} species (conserved)")
    print(f"  - Value = -1 in <= {args.max_missing} species (not too many missing)")
    print()
    
    mask = df.apply(lambda row: filter_row(row, present_species, 
                                            min_species=args.min_conserved,  # Not used but kept for compatibility
                                            min_conserved=args.min_conserved, 
                                            max_missing=args.max_missing,
                                            threshold=args.threshold), axis=1)
    df_filtered = df[mask].copy()
    
    print(f"  Original rows: {len(df)}")
    print(f"  Filtered rows: {len(df_filtered)}")
    print(f"  Rows removed: {len(df) - len(df_filtered)}")
    print(f"  Retention rate: {100 * len(df_filtered) / len(df):.2f}%")
    print()
    
    # Replace species IDs with common names
    print("Replacing species IDs with common names...")
    df_renamed = df_filtered.copy()
    
    # Rename columns
    rename_dict = {}
    for col in df_renamed.columns:
        if col in id_to_name:
            rename_dict[col] = id_to_name[col]
        elif col in TARGET_SPECIES.values():
            # If not in mapping but in target species, use the name from TARGET_SPECIES
            name = [k for k, v in TARGET_SPECIES.items() if v == col][0]
            rename_dict[col] = name
    
    df_renamed = df_renamed.rename(columns=rename_dict)
    
    renamed_count = len([c for c in rename_dict.values() if c in df_renamed.columns])
    print(f"  Renamed {renamed_count} species columns")
    print()
    
    # Format numeric columns to preserve decimal format (no scientific notation)
    print("Formatting numeric columns to preserve decimal format...")
    numeric_cols = [col for col in df_renamed.columns if col != 'decoyID']
    
    def format_number(value):
        """Format number to preserve original format without scientific notation."""
        if pd.isna(value):
            return "-1.000000"
        try:
            float_val = float(value)
            if float_val == -1.0:
                return "-1.000000"
            # Format with exactly 6 decimal places, no scientific notation
            return f"{float_val:.6f}"
        except (ValueError, TypeError):
            return str(value)
    
    for col in numeric_cols:
        df_renamed[col] = df_renamed[col].apply(format_number)
    
    # Write output
    print(f"Writing output to {args.output}...")
    df_renamed.to_csv(args.output, sep='\t', index=False, quoting=3)  # quoting=3 is QUOTE_NONE
    print(f"  Written {len(df_renamed)} rows, {len(df_renamed.columns)} columns")
    print()
    
    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Input rows: {len(df)}")
    print(f"Output rows: {len(df_filtered)}")
    print(f"Filter criteria:")
    print(f"  - Value > {args.threshold} in >= {args.min_conserved} target species (conserved)")
    print(f"  - Value = -1 in <= {args.max_missing} target species (not too many missing)")
    print(f"Target species: {len(target_species_ids)} species")
    print(f"Output file: {args.output}")
    print("=" * 70)

if __name__ == '__main__':
    main()

