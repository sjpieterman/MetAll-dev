#!/bin/bash

# Usage: sheetmaker.sh <directory_to_search>
SEARCH_DIR=$1

if [ -z "$SEARCH_DIR" ]; then
    echo "Usage: $0 <directory_to_search>" >&2
    exit 1
fi

echo "sample,fastq_1,fastq_2,strandedness"

# Pattern 1: Look for *_1_val_1.fq.gz (trimmed paired-end reads)
find "$SEARCH_DIR" -name "*_1_val_1.fq.gz" -o -name "*_1_val_1.fastq.gz" | while read -r R1; do
    R2="${R1/_1_val_1.fq/_2_val_2.fq}"
    R2="${R2/_1_val_1.fastq/_2_val_2.fastq}"

    FILENAME=$(basename "$R1")
    # Extract sample ID (e.g., P23252_102_S1_L003 from P23252_102_S1_L003_1_val_1.fq.gz)
    SAMPLE_ID=$(echo "$FILENAME" | sed 's/_[12]_val_[12]\.f.*$//')

    if [ -f "$R2" ]; then
        ABS_R1=$(readlink -f "$R1")
        ABS_R2=$(readlink -f "$R2")
        echo "${SAMPLE_ID},${ABS_R1},${ABS_R2},unstranded"
    fi
done

# Pattern 2: Look for *_R1_001.fastq.gz (standard Illumina paired-end)
find "$SEARCH_DIR" -name "*_R1_001.fastq.gz" | while read -r R1; do
    R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

    FILENAME=$(basename "$R1")
    BASE_ID=$(echo "$FILENAME" | sed 's/_S[0-9].*//')
    PARENT_FOLDER=$(basename "$(dirname "$R1")")

    if [ "$PARENT_FOLDER" != "$BASE_ID" ] && [ "$PARENT_FOLDER" != "test_data" ] && [ "$PARENT_FOLDER" != "." ]; then
        SAMPLE_ID="${BASE_ID}_${PARENT_FOLDER}"
    else
        SAMPLE_ID=$(echo "$FILENAME" | sed 's/_R1_001.fastq.gz//')
    fi

    if [ -f "$R2" ]; then
        ABS_R1=$(readlink -f "$R1")
        ABS_R2=$(readlink -f "$R2")
        echo "${SAMPLE_ID},${ABS_R1},${ABS_R2},unstranded"
    fi
done