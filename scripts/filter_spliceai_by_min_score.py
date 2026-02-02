#!/usr/bin/env python3
"""
Filter SpliceAI_MAF_protein_coding.tsv to only keep decoyIDs found in DecoySpliceSites_f50.01min.bed
by matching coordinates to decoyIDs.bed (testing both exact match and offset by 24)
"""

import sys

# Read DecoySpliceSites_f50.01min.bed (coordinates with min score > 0.01)
print("Loading DecoySpliceSites_f50.01min.bed...")
decoy_sites_coords = set()
    with open('DecoySpliceSites_f50.01min.bed', 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
            chr_name = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            decoy_sites_coords.add((chr_name, start, end))

print(f"Loaded {len(decoy_sites_coords)} sites from DecoySpliceSites_f50.01min.bed")

# Read decoyIDs.bed (mapping coordinates to decoyIDs)
print("Loading decoyIDs.bed...")
decoy_id_ranges = []  # List of (chr, start, end, decoyID)
decoy_id_ranges_offset = []  # List with offset coordinates for testing
with open('decoyIDs.bed', 'r') as f:
    for line in f:
        fields = line.strip().split()
        if len(fields) >= 4:
            chr_name = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            decoy_id = fields[3]
            decoy_id_ranges.append((chr_name, start, end, decoy_id))
            # Also create offset ranges for testing
            decoy_id_ranges_offset.append((chr_name, start - 24, end - 24, decoy_id))
            decoy_id_ranges_offset.append((chr_name, start + 24, end + 24, decoy_id))

print(f"Loaded {len(decoy_id_ranges)} entries from decoyIDs.bed")

# Test matching: check if DecoySpliceSite coordinates fall within decoyID ranges
matching_decoy_ids_exact = set()
matching_decoy_ids_offset = set()

# Test exact coordinate overlap (DecoySpliceSite within decoyID range)
for site_chr, site_start, site_end in decoy_sites_coords:
    for decoy_chr, decoy_start, decoy_end, decoy_id in decoy_id_ranges:
        if (site_chr == decoy_chr and 
            site_start >= decoy_start and 
            site_end <= decoy_end):
            matching_decoy_ids_exact.add(decoy_id)
            break

# Test offset by 24 coordinate overlap
for site_chr, site_start, site_end in decoy_sites_coords:
    for decoy_chr, decoy_start, decoy_end, decoy_id in decoy_id_ranges_offset:
        if (site_chr == decoy_chr and 
            site_start >= decoy_start and 
            site_end <= decoy_end):
            matching_decoy_ids_offset.add(decoy_id)
            break

print(f"\nMatching results:")
print(f"  Exact overlap: {len(matching_decoy_ids_exact)} decoyIDs")
print(f"  Offset by 24 overlap: {len(matching_decoy_ids_offset)} decoyIDs")

# Use the method that found more matches
if len(matching_decoy_ids_offset) > len(matching_decoy_ids_exact):
    matching_decoy_ids = matching_decoy_ids_offset
    match_method = "offset by 24 overlap"
else:
    matching_decoy_ids = matching_decoy_ids_exact
    match_method = "exact overlap"

print(f"\nUsing {match_method} matching: {len(matching_decoy_ids)} decoyIDs")
    
# Read SpliceAI TSV file and filter
print("\nLoading SpliceAI_MAF_protein_coding.tsv...")
header = None
filtered_lines = []
matched_count = 0
total_count = 0

with open('SpliceAI_MAF_protein_coding.tsv', 'r') as f:
    # Read header
    header = f.readline().strip()
    filtered_lines.append(header)
    
    # Read and filter data
    for line in f:
        total_count += 1
        fields = line.strip().split('\t')
        if len(fields) > 0:
            decoy_id = fields[0]
            if decoy_id in matching_decoy_ids:
                filtered_lines.append(line.strip())
                matched_count += 1

print(f"Loaded {total_count} rows from SpliceAI_MAF_protein_coding.tsv")
print(f"Filtered to {matched_count} rows")

# Save filtered file
    output_file = 'SpliceAI_MAF_protein_coding_min0.01.tsv'
print(f"\nSaving filtered data to {output_file}...")
with open(output_file, 'w') as f:
    f.write('\n'.join(filtered_lines) + '\n')

print(f"Done! Saved {len(filtered_lines) - 1} data rows (plus header)")
