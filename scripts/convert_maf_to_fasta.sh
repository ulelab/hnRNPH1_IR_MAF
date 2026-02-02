#!/bin/bash
# Convert all MAF files in chr2_part5 directory to FASTA format using msaview

INPUT_DIR="chr2_part5"
OUTPUT_DIR="${INPUT_DIR}_fasta"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Count total files
total=$(find "$INPUT_DIR" -name "*.maf" | wc -l)
count=0

echo "Converting MAF files to FASTA format..."
echo "Total files to process: $total"
echo ""

# Process each MAF file
for maf_file in "$INPUT_DIR"/*.maf; do
    if [ -f "$maf_file" ]; then
        count=$((count + 1))
        basename=$(basename "$maf_file" .maf)
        output_file="${OUTPUT_DIR}/${basename}.fa"
        
        echo "[$count/$total] Converting: $(basename "$maf_file")"
        
        # Convert MAF to FASTA using msaview
        msa_view "$maf_file" --in-format MAF --out-format FASTA > "$output_file" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "  ✓ Created: $output_file"
        else
            echo "  ✗ Error converting: $maf_file"
        fi
    fi
done

echo ""
echo "Conversion complete!"
echo "Output directory: $OUTPUT_DIR"
echo "Total files converted: $count"
