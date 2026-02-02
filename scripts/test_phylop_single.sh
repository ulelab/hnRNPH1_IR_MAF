#!/bin/bash
# Test script to process one folder with phyloP

# Paths - adjust these for local testing
MAF_OUTPUT_BASE="/home/mikej10/advbfx/intronretention/maf_output"
OUTPUT_DIR="/home/mikej10/advbfx/intronretention/phylop_test_output"
SCRIPT_DIR="/home/mikej10/advbfx/intronretention"

# Try to find phastCons model file
SCRIPT_DIR_ABS="$(cd "$(dirname "$0")" && pwd)"
if [ -f "${SCRIPT_DIR_ABS}/hg38.phastCons100way.mod" ]; then
    PHYLO_CONS_MOD="${SCRIPT_DIR_ABS}/hg38.phastCons100way.mod"
elif [ -f "/scratch/prj/ppn_rnp_networks/users/mike.jones/data/phast/hg38.phastCons100way.mod" ]; then
    PHYLO_CONS_MOD="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/phast/hg38.phastCons100way.mod"
elif [ -f "${HOME}/data/phast/hg38.phastCons100way.mod" ]; then
    PHYLO_CONS_MOD="${HOME}/data/phast/hg38.phastCons100way.mod"
elif [ -f "./hg38.phastCons100way.mod" ]; then
    PHYLO_CONS_MOD="./hg38.phastCons100way.mod"
else
    echo "ERROR: phastCons model file not found. Please set PHYLO_CONS_MOD environment variable or update the script."
    echo "Expected locations:"
    echo "  - ${SCRIPT_DIR_ABS}/hg38.phastCons100way.mod"
    echo "  - /scratch/prj/ppn_rnp_networks/users/mike.jones/data/phast/hg38.phastCons100way.mod"
    echo "  - ${HOME}/data/phast/hg38.phastCons100way.mod"
    echo "  - ./hg38.phastCons100way.mod"
    exit 1
fi

echo "Using phastCons model: $PHYLO_CONS_MOD"

# Activate conda environment
# Try HPC path first, then local conda
if [ -f "/scratch/prj/ppn_rnp_networks/users/mike.jones/software/mambaforge/etc/profile.d/conda.sh" ]; then
    CONDA_BASE="/scratch/prj/ppn_rnp_networks/users/mike.jones/software/mambaforge"
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate phast-env
elif command -v conda &> /dev/null; then
    # Use local conda if available
    eval "$(conda shell.bash hook)"
    conda activate phast-env 2>/dev/null || echo "WARNING: Could not activate phast-env, trying to use phyloP from PATH"
else
    echo "WARNING: Conda not found, trying to use phyloP from PATH"
fi

# Check if phyloP is available
if ! command -v phyloP &> /dev/null; then
    echo "ERROR: phyloP command not found. Please activate phast-env environment."
    exit 1
fi

echo "Using phyloP: $(which phyloP)"

# Test with chr1_part1 folder
FOLDER_NAME="chr1_part1"
FOLDER_PATH="${MAF_OUTPUT_BASE}/${FOLDER_NAME}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Testing phyloP processing on folder: ${FOLDER_NAME}"
echo "Folder path: $FOLDER_PATH"
echo "Start time: $(date)"
echo "=========================================="

# Check if folder exists
if [ ! -d "$FOLDER_PATH" ]; then
    echo "ERROR: Folder not found: $FOLDER_PATH"
    exit 1
fi

# Count MAF files in this folder
MAF_FILES=($(find "$FOLDER_PATH" -name "*.maf" -o -name "*.maf.gz" | sort))
NUM_MAF=${#MAF_FILES[@]}

if [ "$NUM_MAF" -eq 0 ]; then
    echo "WARNING: No MAF files found in $FOLDER_PATH"
    exit 0
fi

echo "Found $NUM_MAF MAF files"
echo ""

# Create temporary directory for this task
TEMP_DIR="${OUTPUT_DIR}/temp_${FOLDER_NAME}"
mkdir -p "$TEMP_DIR"

# Process first 5 MAF files as a test
PROCESSED=0
FAILED=0
MAX_TEST=5

echo "Processing first $MAX_TEST files as a test..."
echo ""

for MAF_FILE in "${MAF_FILES[@]:0:$MAX_TEST}"; do
    MAF_BASENAME=$(basename "$MAF_FILE" .gz)
    MAF_BASENAME=$(basename "$MAF_BASENAME" .maf)
    
    echo "Processing: $MAF_BASENAME"
    
    # Output file for phyloP SPH results
    PHYLOP_OUTPUT="${TEMP_DIR}/${MAF_BASENAME}.phylop"
    
    # Run phyloP with --wig-scores to get per-site scores
    # Check if file is gzipped
    if [[ "$MAF_FILE" == *.gz ]]; then
        echo "  Running: phyloP --method SPH --mode CONACC --msa-format MAF --wig-scores (gzipped file)"
        phyloP --method SPH --mode CONACC --msa-format MAF --wig-scores "$PHYLO_CONS_MOD" <(zcat "$MAF_FILE") > "$PHYLOP_OUTPUT" 2>&1
    else
        echo "  Running: phyloP --method SPH --mode CONACC --msa-format MAF --wig-scores"
        phyloP --method SPH --mode CONACC --msa-format MAF --wig-scores "$PHYLO_CONS_MOD" "$MAF_FILE" > "$PHYLOP_OUTPUT" 2>&1
    fi
    
    EXIT_CODE=$?
    
    if [ $EXIT_CODE -ne 0 ] || [ ! -s "$PHYLOP_OUTPUT" ]; then
        echo "  WARNING: phyloP failed or produced empty output for $MAF_BASENAME (exit code: $EXIT_CODE)"
        if [ -f "$PHYLOP_OUTPUT" ]; then
            echo "  Output preview:"
            head -10 "$PHYLOP_OUTPUT"
        fi
        rm -f "$PHYLOP_OUTPUT"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    # Show first few lines of output to verify format
    echo "  SUCCESS: Created phyloP output ($(wc -l < "$PHYLOP_OUTPUT") lines)"
    echo "  First 5 lines of output:"
    head -5 "$PHYLOP_OUTPUT" | sed 's/^/    /'
    echo ""
    
    PROCESSED=$((PROCESSED + 1))
done

echo ""
echo "Processed: $PROCESSED files"
echo "Failed: $FAILED files"
echo ""

# Process all phyloP outputs and create summary TSV
if [ $PROCESSED -gt 0 ]; then
    echo "Creating summary TSV..."
    python3 "${SCRIPT_DIR}/parse_phylop_sph.py" "$TEMP_DIR" "${OUTPUT_DIR}/${FOLDER_NAME}_phylop_summary.tsv"
    
    if [ -f "${OUTPUT_DIR}/${FOLDER_NAME}_phylop_summary.tsv" ]; then
        echo ""
        echo "Summary TSV created:"
        head -10 "${OUTPUT_DIR}/${FOLDER_NAME}_phylop_summary.tsv"
    fi
fi

# Keep temp directory for inspection
echo ""
echo "=========================================="
echo "Test completed"
echo "Temporary files: $TEMP_DIR"
echo "Output TSV: ${OUTPUT_DIR}/${FOLDER_NAME}_phylop_summary.tsv"
echo "End time: $(date)"
echo "=========================================="

exit 0

