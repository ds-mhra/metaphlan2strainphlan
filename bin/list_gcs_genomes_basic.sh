#!/bin/bash

# Check if less than 2 arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <bucket_path> <output_csv> . Please note the directory containing NCBI reference files should end with `_ncbi_dataset`"
    exit 1
fi

BUCKET_PATH="$1"
OUTPUT_CSV="$2"

# Create .csv header
echo "ref_genome_strainID,ref_genome_path,species_name,taxonID,associated_clade" > "$OUTPUT_CSV"

# List and process files from GCS
gsutil ls "${BUCKET_PATH}/**/*.fna" | while read -r file; do
    
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
