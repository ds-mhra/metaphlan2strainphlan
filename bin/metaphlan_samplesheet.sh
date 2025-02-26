#!/bin/bash

# Set strict error handling
set -euo pipefail

# Check if less than 2 arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <fastq_dir> <samplesheet_csv> . Ensure input directory contains fastq files that contains `_R1_001` and/or `_R2_001` "
    exit 1
fi

# Define inputs and outdir
FASTQ_DIR="$1"
OUTPUT_CSV="$2"

# Define log file for tracking and clear previous entry
LOG_FILE="samplesheet_creation.log"
> "$LOG_FILE"

# Create samplesheet.csv
# Print header and count samples
echo "sample,fastq_1,fastq_2" > "$OUTPUT_CSV"


# Detect if input is a GCS bucket or local path and list files
if [[ "$FASTQ_DIR" == gs://* ]]; then
    FASTQ_LIST=$(gsutil ls "$FASTQ_DIR" | grep -E "_R[12]_001\.fastq\.gz$")
else
    FASTQ_LIST=$(find "$FASTQ_DIR" -type f -name "*_R[12]_001.fastq.gz")
fi


# Process fastq files
echo "$FASTQ_LIST" | grep "_R1_001\.fastq\.gz$" | while read -r r1; do
    sample=$(basename "$r1" | sed 's/_R1_001\.fastq\.gz//')
    
    # Check if R2 file exists
    r2=$(echo "$r1" | sed 's/_R1_/_R2_/')
    if echo "$FASTQ_LIST" | grep -q "$r2"; then
        echo "${sample},${r1},${r2}" >> "$OUTPUT_CSV"
        echo "Paired sample: $sample" >> "$LOG_FILE"
    else
        # Optional: Add single-end samples
        echo "${sample},${r1}," >> "$OUTPUT_CSV"
        echo "Unpaired sample: $sample" >> "$LOG_FILE"
    fi
done

