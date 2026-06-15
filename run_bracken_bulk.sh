#!/bin/bash
# run_bracken_bulk.sh
# Automates Bracken analysis for Kraken-GPU output reports.

set -euo pipefail

# --- Configuration ---
# Directory containing Kraken-GPU .report.txt files
KRAKEN_DIR="/media/baseripper/volume_28TB_1/Gkraker"

# Path to Bracken executable
BRACKEN_EXE="/home/baseripper/programs/PlayTools/Bracken/bracken"

# Kraken2 Database directory (must contain .kmer_distrib files)
DB_DIR="/media/baseripper/volume_28TB_2/k2_standard_16_GB_20260226"

# Parameters
READ_LEN=150
LEVEL="S"  # S=Species, G=Genus, F=Family, etc.

# Output directory for Bracken results
OUTPUT_DIR="${KRAKEN_DIR}/bracken"
# ---------------------

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

echo "=========================================================="
echo "Starting Bulk Bracken Analysis"
echo "Input:    $KRAKEN_DIR"
echo "Output:   $OUTPUT_DIR"
echo "DB:       $DB_DIR"
echo "Level:    $LEVEL"
echo "Read Len: $READ_LEN"
echo "=========================================================="

# Count reports
REPORT_COUNT=$(ls -1 "$KRAKEN_DIR"/*.report.txt 2>/dev/null | wc -l)
if [ "$REPORT_COUNT" -eq 0 ]; then
    echo "Error: No .report.txt files found in $KRAKEN_DIR"
    exit 1
fi

echo "Found $REPORT_COUNT reports to process."

for report in "$KRAKEN_DIR"/*.report.txt; do
    [ -e "$report" ] || continue
    
    filename=$(basename "$report")
    # Strip .report.txt extension
    sample_name="${filename%.report.txt}"
    
    echo "----------------------------------------------------------"
    echo "Processing Sample: $sample_name"
    
    "$BRACKEN_EXE" \
        -d "$DB_DIR" \
        -i "$report" \
        -o "${OUTPUT_DIR}/${sample_name}_bracken.txt" \
        -r "$READ_LEN" \
        -l "$LEVEL" \
        -w "${OUTPUT_DIR}/${sample_name}_bracken_report.txt"

    echo "Completed $sample_name"
done

echo "=========================================================="
echo "Bracken analysis complete."
echo "Results (standard): ${OUTPUT_DIR}/*_bracken.txt"
echo "Reports (re-estimated): ${OUTPUT_DIR}/*_bracken_report.txt"
echo "=========================================================="
