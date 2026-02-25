#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --array=1-50
#SBATCH --job-name=mafsInRegion_remote
#SBATCH --output=mafsInRegion_remote_%A_%a.out
#SBATCH --error=mafsInRegion_remote_%A_%a.err

# Paths - adjust these for your HPC
STREAM_MAF_SCRIPT="/scratch/prj/ppn_rnp_networks/users/mike.jones/software/stream_maf_regions.py"
BED_PARTS_DIR="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/maf/decoyIDs"
OUTPUT_DIR="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/maf/maf_output_remote"
CHRPARTS_FILE="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/maf/chrparts.txt"
BASE_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/maf"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Get the bed file name for this array task
BED_PART=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHRPARTS_FILE")

if [ -z "$BED_PART" ]; then
    echo "ERROR: No bed file found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

BED_FILE="${BED_PARTS_DIR}/${BED_PART}"

# Check if bed file exists
if [ ! -f "$BED_FILE" ]; then
    echo "ERROR: BED file not found: $BED_FILE"
    exit 1
fi

# Extract chromosome name from bed file name (e.g., chr1_part1.bed -> chr1)
CHR_NAME=$(echo "$BED_PART" | sed 's/_part[0-9]*\.bed$//')

# Output directory for this bed file part (will create separate files per name)
OUTPUT_SUBDIR="${OUTPUT_DIR}/${BED_PART%.bed}"

echo "=========================================="
echo "Task ${SLURM_ARRAY_TASK_ID}: Processing ${BED_PART}"
echo "BED file: $BED_FILE"
echo "Chromosome: $CHR_NAME"
echo "Output directory: $OUTPUT_SUBDIR"
echo "Start time: $(date)"
echo "=========================================="

# Count sites in bed file
NUM_SITES=$(wc -l < "$BED_FILE")
echo "Number of sites: $NUM_SITES"
echo ""

# Create output subdirectory
mkdir -p "$OUTPUT_SUBDIR"

# Stream MAF from remote server and extract regions
python3 "$STREAM_MAF_SCRIPT" \
    "$BED_FILE" \
    --out-dir "$OUTPUT_SUBDIR" \
    --chrom "$CHR_NAME" \
    --base-url "$BASE_URL"

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Task ${SLURM_ARRAY_TASK_ID} completed"
    echo "End time: $(date)"
    
    # Count output files created
    if [ -d "$OUTPUT_SUBDIR" ]; then
        OUTPUT_COUNT=$(find "$OUTPUT_SUBDIR" -name "*.maf" | wc -l)
        OUTPUT_SIZE=$(du -sh "$OUTPUT_SUBDIR" | cut -f1)
        echo "Output files created: $OUTPUT_COUNT"
        echo "Output directory size: $OUTPUT_SIZE"
    fi
else
    echo "ERROR: Task ${SLURM_ARRAY_TASK_ID} failed with exit code $EXIT_CODE"
fi
echo "=========================================="

exit $EXIT_CODE