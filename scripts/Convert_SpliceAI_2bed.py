#!/usr/bin/env python3
"""
Convert SpliceAI output TSV to strand-aware bedgraph format.

Input format: 
  - TSV file with decoyID in first column, then species columns
  - Each species column contains: [max1][pos1][max2][pos2]
    where max1/max2 are splice scores and pos1/pos2 are positions
  
Coordinates:
  - Come from BED file using decoyID match in column 4
  - For positive strand: genomic_pos = start + pos
  - For negative strand: genomic_pos = end - pos
  
Output:
  - bedgraph format (chrom, start, end, score)
  - Scores are signed: positive for + strand, negative for - strand
  - Duplicate positions are averaged
"""

import re
import pandas as pd
import sys

# Input settings - modify these paths as needed
input_path = "/home/mikej10/advbfx/intronretention/all_spliceai_scores_SCN_filtered.tsv"
bed_path = "/home/mikej10/advbfx/intronretention/DecoysInConservedGenes50thpercentile.bed"
output_path = "/home/mikej10/advbfx/intronretention/all_spliceai_scores_SCN_filtered.bedgraph"
logfile = "/home/mikej10/advbfx/intronretention/spliceai_bed_conversion.log"

# Species to process (will look for "Human" or "hg38" in header)
target_species = "Human"  # Options: "Human" or "hg38"

# Logging
log = open(logfile, "w")
bedgraph_records = []

# Load BED file: create mapping from decoyID to (chrom, start, end, strand)
decoy_coords = {}
with open(bed_path, "r") as bed_file:
    for line in bed_file:
        line = line.strip()
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) >= 6:
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            decoy_id = fields[3]
            strand = fields[5]
            decoy_coords[decoy_id] = (chrom, start, end, strand)

print(f"Loaded {len(decoy_coords)} decoy coordinates from BED file")

# Process TSV file
with open(input_path, "r") as infile:
    header_line = infile.readline().strip()
    header_fields = header_line.split("\t")
    
    if len(header_fields) < 2:
        log.write(f"ERROR: Invalid header line: {header_line}\n")
        sys.exit(1)
    
    # First column is decoyID, rest are species names
    species_names = header_fields[1:]
    print(f"Found {len(species_names)} species columns")
    
    # Find the target species column (Human or hg38)
    target_species_idx = None
    for idx, species_name in enumerate(species_names):
        if species_name == target_species or species_name == "hg38":
            target_species_idx = idx
            print(f"Found target species column: '{species_name}' at index {idx}")
            break
    
    if target_species_idx is None:
        log.write(f"ERROR: Could not find '{target_species}' or 'hg38' in header columns\n")
        print(f"ERROR: Could not find '{target_species}' or 'hg38' in header columns")
        print(f"Available columns: {species_names[:10]}...")
        sys.exit(1)
    
    # Process each data line
    for line_num, line in enumerate(infile, start=2):
        line = line.strip()
        if not line:
            continue
        
        fields = line.split("\t")
        if len(fields) < target_species_idx + 2:  # Need decoyID + target species column
            log.write(f"WARNING: Skipping line {line_num}: insufficient fields\n")
            continue
        
        decoy_id = fields[0]
        
        # Look up coordinates for this decoy
        if decoy_id not in decoy_coords:
            log.write(f"WARNING: Decoy ID not found in BED file (line {line_num}): {decoy_id}\n")
            continue
        
        chrom, start_coord, end_coord, strand = decoy_coords[decoy_id]
        
        # Process only the target species column
        species_data = fields[target_species_idx + 1]
        
        # Parse [max1][pos1][max2][pos2] format
        # Pattern: [float][int][float][int]
        match = re.match(r'\[([\d.]+)\]\[(\d+)\]\[([\d.]+)\]\[(\d+)\]', species_data)
        if not match:
            # Check if it's a valid format (might have -1.0 for missing data)
            if '[' in species_data and ']' in species_data:
                log.write(f"WARNING: Could not parse species data for {decoy_id}, {target_species} (line {line_num}): {species_data}\n")
            continue
        
        max1 = float(match.group(1))
        pos1 = int(match.group(2))
        max2 = float(match.group(3))
        pos2 = int(match.group(4))
        
        # Skip if scores are -1.0 (missing/invalid data)
        if max1 == -1.0 and max2 == -1.0:
            continue
        
        # Convert positions to genomic coordinates based on strand
        if strand == "+":
            # Position is nucleotides from start coordinate
            genomic_pos1 = start_coord + pos1
            genomic_pos2 = start_coord + pos2
            signed_score1 = max1
            signed_score2 = max2
        else:  # strand == "-"
            # Position is nucleotides BEFORE end coordinate
            genomic_pos1 = end_coord - pos1
            genomic_pos2 = end_coord - pos2
            # Negative strand: make scores negative
            signed_score1 = -max1
            signed_score2 = -max2
        
        # Add bedgraph entries (one for each score/position pair)
        # Only add if score is valid (> -1.0)
        if max1 != -1.0:
            bedgraph_records.append([chrom, genomic_pos1, genomic_pos1 + 1, signed_score1])
        if max2 != -1.0:
            bedgraph_records.append([chrom, genomic_pos2, genomic_pos2 + 1, signed_score2])

log.close()

# Convert to pandas for sorting and collapsing
if len(bedgraph_records) == 0:
    print("WARNING: No bedgraph records generated!")
    sys.exit(1)

df = pd.DataFrame(bedgraph_records, columns=["chrom", "start", "end", "score"])
df = df.sort_values(by=["chrom", "start", "end"])

# Collapse duplicates by averaging signal (if same position appears multiple times)
collapsed = df.groupby(["chrom", "start", "end"], as_index=False).agg({"score": "mean"})

# Output final bedGraph
collapsed.to_csv(output_path, sep='\t', header=False, index=False, float_format='%.6f')

print(f"✅ Strand-aware bedGraph written: {output_path}")
print(f"   Total records: {len(collapsed)}")
print(f"   Log file: {logfile}")
