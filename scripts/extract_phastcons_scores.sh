#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Usage:
#   bash extract_phastcons_scores.sh <bigwig_file> <bed_file>
# Output (stdout):
#   name<TAB>mean_phastCons

BIGWIG_FILE="${1:?Usage: bash extract_phastcons_scores.sh <bigwig_file> <bed_file>}"
BED_FILE="${2:?Usage: bash extract_phastcons_scores.sh <bigwig_file> <bed_file>}"

# Resolve bigWigAverageOverBed binary without hardcoded absolute path.
if command -v bigWigAverageOverBed >/dev/null 2>&1; then
    BINARY_PATH="$(command -v bigWigAverageOverBed)"
elif [ -x "$PROJECT_DIR/bin/bigWigAverageOverBed" ]; then
    BINARY_PATH="$PROJECT_DIR/bin/bigWigAverageOverBed"
elif [ -x "$PROJECT_DIR/bin/$(uname -m)/bigWigAverageOverBed" ]; then
    BINARY_PATH="$PROJECT_DIR/bin/$(uname -m)/bigWigAverageOverBed"
elif [ -x "$HOME/bin/$(uname -m)/bigWigAverageOverBed" ]; then
    BINARY_PATH="$HOME/bin/$(uname -m)/bigWigAverageOverBed"
else
    echo "Error: bigWigAverageOverBed not found in PATH or project bin locations." >&2
    exit 1
fi

if [ ! -f "$BIGWIG_FILE" ]; then
    echo "Error: BigWig file not found: $BIGWIG_FILE" >&2
    exit 1
fi

if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file not found: $BED_FILE" >&2
    exit 1
fi

TEMP_BED="$(mktemp)"
TEMP_OUT="$(mktemp)"

# Keep only chr/start/end/name and enforce unique names.
awk -F'\t' '
BEGIN {OFS="\t"}
{
    name = $4
    if (!seen[name]) {
        print $1, $2, $3, name
        seen[name] = 1
    }
}
' "$BED_FILE" > "$TEMP_BED"

"$BINARY_PATH" "$BIGWIG_FILE" "$TEMP_BED" "$TEMP_OUT"

# bigWigAverageOverBed output columns:
# name, size, covered, sum, mean0, mean
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $6}' "$TEMP_OUT"

rm -f "$TEMP_BED" "$TEMP_OUT"
