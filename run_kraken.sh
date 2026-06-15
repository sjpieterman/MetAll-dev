#!/bin/bash
set -euo pipefail

# Add kraken2 to PATH if not already there
export PATH="/home/baseripper/programs/kraken2:$PATH"

input_folder="/media/baseripper/volume_28TB_1/01_RawDataFull_P11065/results_run_01_NoBracken/unmapped_reads/unmapped"
output_folder="/media/baseripper/volume_28TB_1/01_RawDataFull_P11065/results_run_01_NoBracken/kraken_vanlig"
kraken2_db="/media/baseripper/volume_28TB_2/k2_standard_16_GB_20260226"
threads=56

mkdir -p "${input_folder}/kraken"
mkdir -p "${input_folder}/bracken"
mkdir -p "${output_folder}/kraken"
mkdir -p "${output_folder}/bracken"

# Find paired uncompressed FASTQ files in the test_files folder.
find "$input_folder" -type f -name "*.mate1.gz" | while IFS= read -r r1; do
    r2="${r1%.mate1.gz}.mate2.gz"

    if [[ ! -f "$r2" ]]; then
        echo "WARNING: missing R2 for: $r1 (expected: $r2). Skipping."
        continue
    fi

    filename="$(basename "$r1")"
    basename="${filename%.mate1.fastq}"

    echo "Running Kraken2 for $basename..."

    kraken2 --db "$kraken2_db" --threads "$threads" \
        --report "${input_folder}/kraken/${basename}_kraken2.tax" \
        --paired "$r1" "$r2" \
        --output "${output_folder}/kraken/${basename}_kraken2.krk"

    #/home/baseripper/programs/Bracken/bracken \
    #    -d "$kraken2_db" \
    #    -i "${input_folder}/kraken/${basename}_kraken2.tax" \
    #    -o "${output_folder}/bracken/${basename}_bracken.txt" \
    #    -r 150 -l S
done

echo "Processing complete. Kraken2 and Bracken analyses have been saved."
