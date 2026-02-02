#!/usr/bin/env python3
"""
Find top 50 max scores from SpliceAI_MAF_non_exonic_no_unstranded_overlaps.tsv and convert to BED format.
"""

import pandas as pd
import sys

def main():
    # File paths
    spliceai_file = "SpliceAI_MAF_non_exonic_no_unstranded_overlaps.tsv"
    decoy_bed_file = "decoyIDs.bed"
    output_bed = "SpliceAI_MAF_non_exonic_no_unstranded_overlaps_top50.bed"
    
    print("=" * 70)
    print("Extracting Top 50 Max SpliceAI Scores (No Unstranded Overlaps)")
    print("=" * 70)
    print()
    
    print("Step 1: Loading SpliceAI MAF data...")
    # Read SpliceAI data as strings to preserve format
    spliceai_df = pd.read_csv(spliceai_file, sep="\t", dtype=str)
    print(f"  Loaded {len(spliceai_df)} rows")
    
    # Get all species columns (all columns except decoyID)
    species_cols = [col for col in spliceai_df.columns if col != 'decoyID']
    print(f"  Found {len(species_cols)} species columns")
    print()
    
    print("Step 2: Calculating max scores for each decoy ID...")
    # Convert to numeric for calculation, handling -1.000000 as NaN
    # Keep decoyID as string
    numeric_df = spliceai_df[['decoyID']].copy()
    
    # Convert species columns to numeric
    for col in species_cols:
        numeric_df[col] = pd.to_numeric(spliceai_df[col], errors='coerce')
        # Replace -1.0 with NaN
        numeric_df[col] = numeric_df[col].replace(-1.0, pd.NA)
    
    # Calculate max score across all species for each row
    numeric_df['max_score'] = numeric_df[species_cols].max(axis=1, skipna=True)
    
    # Remove rows where max_score is NaN (all -1 values)
    numeric_df = numeric_df.dropna(subset=['max_score'])
    print(f"  Rows with valid scores: {len(numeric_df)}")
    print()
    
    # Sort by max_score descending and take top 50
    top50_df = numeric_df.sort_values('max_score', ascending=False).head(50)[['decoyID', 'max_score']].copy()
    print("Step 3: Top 50 max scores:")
    print(f"  Min max_score: {top50_df['max_score'].min():.6f}")
    print(f"  Max max_score: {top50_df['max_score'].max():.6f}")
    print(f"  Mean max_score: {top50_df['max_score'].mean():.6f}")
    print()
    
    print(f"Step 4: Loading decoy BED file: {decoy_bed_file}")
    # Load decoy BED file to map decoyID to coordinates
    # Try tab-separated first, then space-separated
    try:
        decoy_bed_df = pd.read_csv(decoy_bed_file, sep='\t', header=None,
                                   names=['chr', 'start', 'end', 'name', 'score', 'strand', 'extra'],
                                   usecols=[0, 1, 2, 3, 4, 5])
        print(f"  Loaded {len(decoy_bed_df)} decoy entries (tab-separated)")
    except:
        try:
            decoy_bed_df = pd.read_csv(decoy_bed_file, sep=r'\s+', header=None,
                                       names=['chr', 'start', 'end', 'name', 'score', 'strand', 'extra'],
                                       usecols=[0, 1, 2, 3, 4, 5])
            print(f"  Loaded {len(decoy_bed_df)} decoy entries (space-separated)")
        except Exception as e:
            print(f"ERROR: Failed to read decoy BED file: {e}")
            sys.exit(1)
    print()
    
    print(f"Step 5: Matching decoy IDs and creating BED file...")
    # Merge top50 with decoy BED to get coordinates
    merged_df = top50_df.merge(decoy_bed_df, left_on='decoyID', right_on='name', how='inner')
    
    if len(merged_df) < len(top50_df):
        missing = len(top50_df) - len(merged_df)
        print(f"  WARNING: {missing} decoy IDs not found in decoy BED file")
    
    print(f"  Matched {len(merged_df)} entries")
    print()
    
    # Create BED output with max_score as the score column
    bed_output = merged_df[['chr', 'start', 'end', 'decoyID', 'max_score', 'strand']].copy()
    
    # Format max_score to 6 decimal places
    bed_output['max_score'] = bed_output['max_score'].apply(lambda x: f"{x:.6f}")
    
    # Sort by max_score descending
    bed_output = bed_output.sort_values('max_score', ascending=False, key=lambda x: pd.to_numeric(x))
    
    print(f"Step 6: Writing BED file: {output_bed}")
    # Write BED file (tab-separated, no header)
    bed_output.to_csv(output_bed, sep='\t', index=False, header=False, quoting=3)
    
    print(f"  Written {len(bed_output)} entries")
    print()
    print("Top 5 entries:")
    print(bed_output.head().to_string(index=False))
    print()
    
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Input file: {spliceai_file}")
    print(f"Total decoy IDs: {len(spliceai_df)}")
    print(f"Top 50 max scores extracted")
    print(f"Output BED file: {output_bed}")
    print(f"BED entries: {len(bed_output)}")
    print("=" * 70)

if __name__ == "__main__":
    main()
