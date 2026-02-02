#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=1-53
#SBATCH --job-name=phylop_array
#SBATCH --output=phylop_%A_%a.out
#SBATCH --error=phylop_%A_%a.err

# Paths
MAF_OUTPUT_BASE="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/maf/maf_output"
PHYLO_CONS_MOD="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/phast/hg38.phastCons100way.mod"
OUTPUT_DIR="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/phylop_output"
SCRIPT_DIR="/scratch/prj/ppn_rnp_networks/users/mike.jones/scripts"

# Activate conda environment
CONDA_BASE="/scratch/prj/ppn_rnp_networks/users/mike.jones/software/mambaforge"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate phast-env

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Get the folder for this array task
# List all directories in maf_output and select by array index
# Use find to explicitly list only directories (one level deep, excluding base dir)
FOLDERS=($(find "${MAF_OUTPUT_BASE}" -mindepth 1 -maxdepth 1 -type d | sort))
TOTAL_FOLDERS=${#FOLDERS[@]}

if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_FOLDERS" ]; then
    echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} exceeds number of folders (${TOTAL_FOLDERS})"
    exit 1
fi

FOLDER_PATH="${FOLDERS[$((SLURM_ARRAY_TASK_ID - 1))]}"
FOLDER_NAME=$(basename "$FOLDER_PATH")

echo "=========================================="
echo "Task ${SLURM_ARRAY_TASK_ID}: Processing folder ${FOLDER_NAME}"
echo "Folder path: $FOLDER_PATH"
echo "Start time: $(date)"
echo "=========================================="

# Count MAF files in this folder (only .maf files, no gzipped)
MAF_FILES=($(find "$FOLDER_PATH" -name "*.maf" -type f | sort))
NUM_MAF=${#MAF_FILES[@]}

if [ "$NUM_MAF" -eq 0 ]; then
    echo "WARNING: No MAF files found in $FOLDER_PATH"
    exit 0
fi

echo "Found $NUM_MAF MAF files"
echo ""

# Create temporary directory for this task
TEMP_DIR="${OUTPUT_DIR}/temp_${FOLDER_NAME}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TEMP_DIR"

# Process each MAF file
PROCESSED=0
FAILED=0

for MAF_FILE in "${MAF_FILES[@]}"; do
    MAF_BASENAME=$(basename "$MAF_FILE" .maf)
    
    echo "Processing: $MAF_BASENAME"
    
    # Output file for phyloP SPH results
    PHYLOP_OUTPUT="${TEMP_DIR}/${MAF_BASENAME}.phylop"
    
    # Run phyloP with --wig-scores to get per-site scores
    phyloP --method SPH --mode CONACC --msa-format MAF --wig-scores "$PHYLO_CONS_MOD" "$MAF_FILE" > "$PHYLOP_OUTPUT" 2>&1
    
    EXIT_CODE=$?
    
    if [ $EXIT_CODE -ne 0 ] || [ ! -s "$PHYLOP_OUTPUT" ]; then
        echo "  WARNING: phyloP failed or produced empty output for $MAF_BASENAME"
        rm -f "$PHYLOP_OUTPUT"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    PROCESSED=$((PROCESSED + 1))
    
    if [ $((PROCESSED % 100)) -eq 0 ]; then
        echo "  Progress: $PROCESSED/$NUM_MAF files processed"
    fi
done

echo ""
echo "Processed: $PROCESSED files"
echo "Failed: $FAILED files"
echo ""

# Process all phyloP outputs and create summary TSV
echo "Creating summary TSV..."
python3 "${SCRIPT_DIR}/parse_phylop_sph.py" "$TEMP_DIR" "${OUTPUT_DIR}/${FOLDER_NAME}_phylop_summary.tsv"

# Clean up temporary directory
rm -rf "$TEMP_DIR"

echo ""
echo "=========================================="
echo "SUCCESS: Task ${SLURM_ARRAY_TASK_ID} completed"
echo "Output: ${OUTPUT_DIR}/${FOLDER_NAME}_phylop_summary.tsv"
echo "End time: $(date)"
echo "=========================================="

exit 0

