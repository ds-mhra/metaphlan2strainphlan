#!/bin/bash

# Check if less than 2 arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <bucket_path1> [bucket_path2 ...] <output_csv> . Please note the directory containing NCBI reference files should end with `_ncbi_dataset`"
    exit 1
fi

# Extract last argument as output CSV
OUTPUT_CSV="${!#}"
# Extract all arguments except the last one as bucket paths
INPUT_PATHS=("${@:1:$#-1}")

# Create .csv header
echo "ref_genome_strainID,ref_genome_path,species_name,taxonID,associated_clade" > "$OUTPUT_CSV"

# Process each bucket
for INPUT_PATH in "${INPUT_PATHS[@]}"; do

    # Detect if input is a GCS bucket or local path, then list
    if [[ "$INPUT_PATH" == gs://* ]]; then
        FILE_LIST=$(gsutil ls "${INPUT_PATH}/**/*.fna")
    else
        FILE_LIST=$(find "$INPUT_PATH" -type f -name "*.fna")
    fi

    # Process files from GCS
    echo "$FILE_LIST" | while read -r file; do
        
        strain_id=$(basename "$(dirname "$file")")      # OR `| grep -oP 'GCA_\d+\.\d+|GCF_\d+\.\d+')`

        # Try multiple strategies to extract species_name
        species_name=$(echo "$file" | grep -Eo '/[0-9]+_ncbi_dataset/' | tr -d '/_ncbi_dataset' || \
                    echo "$file" | grep -Eo '/[0-9]+/' | tr -d '/')
        species_name=${species_name:-$(basename "$(dirname "$(dirname "$(dirname "$file")")")" | sed "s/_ncbi_dataset//g" | tr '_' ' ')}

        # Debug message
        if [ -z "$species_name" ]; then
            echo "Warning: Could not extract species name for $file" >&2
        fi

        # Write to .csv
        echo "${strain_id},${file},${species_name},," >> "$OUTPUT_CSV"

    done
done
