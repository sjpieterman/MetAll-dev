#!/bin/bash
# run_folder_gpu.sh
# Usage: ./run_folder_gpu.sh

set -euo pipefail

# Configuration
INPUT_DIR="/media/baseripper/volume_28TB_1/01_RawDataFull_P11065/results_run_01_NoBracken/unmapped_reads/unmapped"
OUTPUT_DIR="/media/baseripper/volume_28TB_1/01_RawDataFull_P11065/results_run_01_NoBracken/Gkraker_2"
DB_DIR="/media/baseripper/volume_28TB_2/k2_standard_16_GB_20260226"
BINARY="/media/baseripper/KINGSTON/projects/GKraker/cmake-build-debug/kraken_gpu"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

echo "Starting bulk Kraken-GPU analysis (Multi-threaded I/O)..."
echo "Input:  $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "DB:     $DB_DIR"

# Collect all R1/R2 pairs
FILES=()
while read -r R1_GZ; do
    R2_GZ=$(echo "$R1_GZ" | sed 's/mate1/mate2/')
    if [ -f "$R2_GZ" ]; then
        FILES+=("$R1_GZ" "$R2_GZ")
    else
        echo "WARNING: No R2 file found for: $(basename "$R1_GZ"). Skipping."
    fi
done < <(find "$INPUT_DIR" -name "*.mate1.gz" | sort)

if [ ${#FILES[@]} -eq 0 ]; then
    echo "No valid pairs found. Exiting."
    exit 1
fi

echo "Found ${#FILES[@]} files ( $(( ${#FILES[@]} / 2 )) pairs )."

# Execute the GPU classifier directly with all files.
# Passing all files at once avoids reloading the database for each sample.
"$BINARY" "${FILES[@]}" \
    -o "$OUTPUT_DIR" \
    -d "$DB_DIR" \
    --paired \

echo "=========================================================="
echo "All samples in $INPUT_DIR have been processed."
echo "Results are in $OUTPUT_DIR"
echo "=========================================================="
