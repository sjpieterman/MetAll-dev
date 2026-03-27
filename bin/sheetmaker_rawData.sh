#!/bin/bash

# Usage: sheetmaker.sh <directory_to_search>
SEARCH_DIR=$1

if [ -z "$SEARCH_DIR" ]; then
    echo "Usage: $0 <directory_to_search>" >&2
    exit 1
fi

echo "sample,fastq_1,fastq_2,strandedness"

find "$SEARCH_DIR" -name "*_R1_001.fastq.gz" | while read -r R1; do
    R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

    FILENAME=$(basename "$R1")
    # Base ID (e.g., P29453_456)
    BASE_ID=$(echo "$FILENAME" | sed 's/_S[0-9].*//')

    # Get the name of the folder containing the file
    PARENT_FOLDER=$(basename "$(dirname "$R1")")

    # Logic: If the parent folder has a unique name (like a run ID)
    # that is NOT the same as the base ID, append it to make the sample name unique.
    if [ "$PARENT_FOLDER" != "$BASE_ID" ] && [ "$PARENT_FOLDER" != "test_data" ] && [ "$PARENT_FOLDER" != "." ]; then
        SAMPLE_ID="${BASE_ID}_${PARENT_FOLDER}"
    else
        # Otherwise, use the full filename prefix as requested
        SAMPLE_ID=$(echo "$FILENAME" | sed 's/_R1_001.fastq.gz//')
    fi

    if [ -f "$R2" ]; then
        ABS_R1=$(readlink -f "$R1")
        ABS_R2=$(readlink -f "$R2")
        echo "${SAMPLE_ID},${ABS_R1},${ABS_R2},unstranded"
    fi
done