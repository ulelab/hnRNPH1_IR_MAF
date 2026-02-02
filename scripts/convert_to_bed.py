#!/usr/bin/env python3
"""
Convert filtered SpliceAI MAF data to BED format for IGV visualization.
Uses Human value as the score and looks up chromosome coordinates from decoyIDs.bed.
"""

import pandas as pd
import argparse
import sys

def load_decoy_bed(decoy_bed_file):
    """
    Load decoy ID to chromosome mapping from BED file.
    Returns a dictionary: decoy_id -> (chr, start, end, strand)
    """
    decoy_map = {}
    
    try:
        # Try tab-delimited first
        bed_df = pd.read_csv(decoy_bed_file, sep='\t', header=None,
                            names=['chr', 'start', 'end', 'name', 'score', 'strand'],
                            usecols=[0, 1, 2, 3, 4, 5])
    except:
        # Try space-delimited
        bed_df = pd.read_csv(decoy_bed_file, sep=r'\s+', header=None,
                            names=['chr', 'start', 'end', 'name', 'score', 'strand'],
                            usecols=[0, 1, 2, 3, 4, 5])
    
    for _, row in bed_df.iterrows():
        decoy_id = str(row['name'])
        decoy_map[decoy_id] = {
            'chr': str(row['chr']),
            'start': int(row['start']),
            'end': int(row['end']),
            'strand': str(row['strand'])
        }
    
    return decoy_map

def parse_decoy_id(decoy_id):
    """
    Parse decoy ID like 'AC003958.2_41411483' or 'ACTN4_38711360'.
    Returns (gene_name, coordinate).
    """
    parts = decoy_id.rsplit('_', 1)
    if len(parts) != 2:
        return None, None
    gene_name = parts[0]
    try:
        coordinate = int(parts[1])
        return gene_name, coordinate
    except ValueError:
        return None, None

def main():
    parser = argparse.ArgumentParser(
        description='Convert filtered SpliceAI MAF data to BED format for IGV'
    )
    parser.add_argument('--input', required=True,
                       help='Input filtered TSV file (SpliceAI_MAF_filtered_conservation.tsv)')
    parser.add_argument('--decoy-bed', required=True,
                       help='decoyIDs.bed file with chromosome coordinates')
    parser.add_argument('--output', required=True,
                       help='Output BED file')
    parser.add_argument('--window', type=int, default=50,
                       help='Window size around coordinate for BED (default: 50 bp)')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("Converting to BED Format for IGV")
    print("=" * 70)
    print(f"Input file: {args.input}")
    print(f"Decoy BED file: {args.decoy_bed}")
    print(f"Output file: {args.output}")
    print(f"Window size: ±{args.window} bp")
    print("=" * 70)
    print()
    
    # Load decoy ID mapping
    print("Loading decoy ID to chromosome mapping...")
    decoy_map = load_decoy_bed(args.decoy_bed)
    print(f"  Loaded {len(decoy_map)} decoy IDs")
    print()
    
    # Load filtered data
    print("Loading filtered SpliceAI data...")
    df = pd.read_csv(args.input, sep='\t')
    print(f"  Loaded {len(df)} rows")
    print()
    
    # Check if Human column exists
    if 'Human' not in df.columns:
        print("ERROR: 'Human' column not found in input file", file=sys.stderr)
        sys.exit(1)
    
    # Create BED entries
    print("Creating BED entries...")
    bed_entries = []
    missing_coords = 0
    found_coords = 0
    
    for _, row in df.iterrows():
        decoy_id = str(row['decoyID'])
        human_value = row['Human']
        
        # Format human value to ensure proper decimal format (no scientific notation)
        try:
            if pd.isna(human_value) or human_value == "-1.000000" or human_value == "-1.0" or human_value == "-1":
                human_score = "0.000000"
            else:
                float_val = float(human_value)
                if float_val == -1.0:
                    human_score = "0.000000"
                else:
                    # Format with 6 decimal places, no scientific notation
                    human_score = f"{float_val:.6f}"
        except (ValueError, TypeError):
            human_score = "0.000000"
        
        # Look up chromosome coordinates from decoy BED
        if decoy_id in decoy_map:
            coord_info = decoy_map[decoy_id]
            chr_name = coord_info['chr']
            start = coord_info['start']
            end = coord_info['end']
            strand = coord_info['strand']
            
            # Expand window around the region
            bed_start = max(0, start - args.window)
            bed_end = end + args.window
            
            bed_entries.append({
                'chr': chr_name,
                'start': bed_start,
                'end': bed_end,
                'name': decoy_id,
                'score': human_score,
                'strand': strand
            })
            found_coords += 1
        else:
            # Try to parse from decoy ID if not in mapping
            gene_name, coordinate = parse_decoy_id(decoy_id)
            if coordinate is not None:
                # We have coordinate but need chromosome - this is a problem
                # For now, skip or use a default
                print(f"  WARNING: Decoy ID {decoy_id} not found in decoy BED, coordinate={coordinate} but no chromosome info", file=sys.stderr)
                missing_coords += 1
            else:
                print(f"  WARNING: Could not parse decoy ID: {decoy_id}", file=sys.stderr)
                missing_coords += 1
    
    print(f"  Found coordinates: {found_coords}")
    print(f"  Missing coordinates: {missing_coords}")
    print()
    
    if len(bed_entries) == 0:
        print("ERROR: No BED entries created", file=sys.stderr)
        sys.exit(1)
    
    # Create BED DataFrame
    bed_df = pd.DataFrame(bed_entries)
    
    # Sort by chromosome and start position
    bed_df = bed_df.sort_values(['chr', 'start', 'end'])
    
    # Write BED file
    print(f"Writing BED file to {args.output}...")
    bed_df.to_csv(args.output, sep='\t', header=False, index=False,
                  columns=['chr', 'start', 'end', 'name', 'score', 'strand'])
    print(f"  Written {len(bed_df)} entries")
    print()
    
    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Input rows: {len(df)}")
    print(f"BED entries created: {len(bed_df)}")
    print(f"Missing coordinates: {missing_coords}")
    print(f"Output file: {args.output}")
    print("=" * 70)

if __name__ == '__main__':
    main()

