#!/bin/bash
# Script to add unique IDs to all split BED files
# Usage: ./add_unique_ids_to_all.sh <input_splitbed_dir> <output_splitbed_dir> [chrparts.txt]

INPUT_DIR="${1:-splitbed}"
OUTPUT_DIR="${2:-splitbed_with_ids}"
CHRPARTS_FILE="${3:-chrparts.txt}"

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Adding unique IDs to all BED files"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# If chrparts.txt exists, use it to process files in order
if [ -f "$CHRPARTS_FILE" ]; then
    echo "Processing files from $CHRPARTS_FILE"
    TOTAL_FILES=$(wc -l < "$CHRPARTS_FILE")
    echo "Total files: $TOTAL_FILES"
    echo ""
    
    while IFS= read -r bed_file; do
        INPUT_FILE="${INPUT_DIR}/${bed_file}"
        OUTPUT_FILE="${OUTPUT_DIR}/${bed_file}"
        
        if [ ! -f "$INPUT_FILE" ]; then
            echo "WARNING: File not found: $INPUT_FILE"
            continue
        fi
        
        echo "Processing: $bed_file"
        
        # Add unique IDs: gene_name_original_coordinate (where original = col2 + 24)
        awk -F'\t' '{
            gene_name = $4
            original_start = $2 + 24
            unique_id = gene_name "_" original_start
            $4 = unique_id
            print
        }' "$INPUT_FILE" > "$OUTPUT_FILE"
        
        INPUT_LINES=$(wc -l < "$INPUT_FILE")
        OUTPUT_LINES=$(wc -l < "$OUTPUT_FILE")
        echo "  $INPUT_LINES lines -> $OUTPUT_LINES lines"
        
    done < "$CHRPARTS_FILE"
else
    echo "chrparts.txt not found. Processing all .bed files in $INPUT_DIR"
    
    for INPUT_FILE in "${INPUT_DIR}"/*.bed; do
        if [ ! -f "$INPUT_FILE" ]; then
            continue
        fi
        
        BED_FILE=$(basename "$INPUT_FILE")
        OUTPUT_FILE="${OUTPUT_DIR}/${BED_FILE}"
        
        echo "Processing: $BED_FILE"
        
        # Add unique IDs
        awk -F'\t' '{
            gene_name = $4
            original_start = $2 + 24
            unique_id = gene_name "_" original_start
            $4 = unique_id
            print
        }' "$INPUT_FILE" > "$OUTPUT_FILE"
        
        INPUT_LINES=$(wc -l < "$INPUT_FILE")
        OUTPUT_LINES=$(wc -l < "$OUTPUT_FILE")
        echo "  $INPUT_LINES lines -> $OUTPUT_LINES lines"
    done
fi

echo ""
echo "Done! All files processed."
echo "Output directory: $OUTPUT_DIR"

