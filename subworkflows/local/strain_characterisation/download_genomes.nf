
process DOWNLOAD_GENOMES {
    publishDir "${params.outdir}/genome_downloads", mode: params.publish_dir_mode
    // errorStrategy 'ignore'

    input:
    val taxon

    output:
    tuple val(meta), path ("${taxon}_ncbi_dataset")    , emit: downloaded_genomes_success  , optional: true
    tuple val(meta), path ("no_ncbi_records.txt")      , emit: downloaded_genomes_failed   , optional: true
    tuple val(meta), path ("successful_taxons.txt")    , emit: downloads_list              , optional: true
    tuple val(meta), path ("alt_ncbi_names.txt")       , emit: alt_genome_names            , optional: true
    

    script:
    """
    #!/bin/bash
    set -e

    echo "${taxon} being processed..."

    datasets download genome taxon "${taxon}" --assembly-level complete --dehydrated --filename "${taxon}_ncbi_dataset.zip"; then
    if [ \$? -eq 0 ]; then
        mkdir -p "tmp/${taxon}"
        unzip -n -d "tmp/${taxon}" "${taxon}_ncbi_dataset.zip"
        datasets rehydrate --directory "tmp/"
        rm "${taxon}_ncbi_dataset.zip"
        echo "${taxon}" >> "successful_taxons.txt"
                    #else
                    #    echo "${taxon}" >> "no_ncbi_records.txt"
    fi

    command=\$(datasets download genome taxon "${taxon}" --assembly-level complete --dehydrated --filename "${taxon}_ncbi_dataset.zip" 2>&1)
    [[ ! -f "alt_ncbi_names.txt" ]] && touch alt_ncbi_names.txt
    [[ ! -f "no_ncbi_records.txt" ]] && touch no_ncbi_records.txt
    
    # If the below error message appears, extract the tool's suggested names (and taxids) to alt_ncbi_names.txt
    
    if [[ "\${command}" == *"Try using one of the suggested taxids:"* ]]; then
        echo "\${command}" | awk '/Try using one of the suggested taxids:/, /^\$/' <<< "\${command}" >> "alt_ncbi_names.txt"
    # Redirect taxons w no further suggestions (from ncbi's tool) to no_records.txt (including species that had issues downloading)
    elif [[ "\${command}" == *"Error:"* || "\${command}" == *"i/o timeout"* ]]; then
        echo "${taxon}" >> "no_ncbi_records.txt"
    fi


    """
}


process RETRY_DOWNLOADS {
    publishDir "${params.outdir}/genome_downloads", mode: params.publish_dir_mode
    // errorStrategy 'ignore'

    input:
    val taxon

    output:
    path "${taxon}_ncbi_dataset", emit: success, optional: true
    path "no_records.txt", emit: failed, optional: true

    script:
    """
    #!/bin/bash

    echo "${taxon} being processed for a final time..."

    datasets download genome taxon "${taxon}" --assembly-level complete --dehydrated --filename "${taxon}_ncbi_dataset.zip"; then
    [[ ! -f "successful_taxons.txt" ]] && touch "successful_taxons.txt"
    [[ ! -f "no_records.txt" ]] && touch "no_records.txt"

    if [ \$? -eq 0 ]; then
        mkdir -p "tmp/${taxon}"
        unzip -n -d "tmp/${taxon}" "${taxon}_ncbi_dataset.zip"
        datasets rehydrate --directory "tmp/"
        rm "${taxon}_ncbi_dataset.zip"
        printf '%s\n' "The taxon ID, ${taxon}, was successful. \n"
        echo "${taxon}" >> "successful_taxons.txt"
    else
        printf '%s\n' "The taxon ID, ${taxon}, was NOT successful. \n"
        echo "${taxon}" >> "no_records.txt"
    fi



    #command=\$(datasets download genome taxon "$taxonID" --assembly-level complete --dehydrated --filename "${tempDir}/${taxonID}_ncbi_dataset.zip" 2>&1)
    #datasets download genome taxon "$taxonID" --assembly-level complete --dehydrated --filename "${tempDir}/${taxonID}_ncbi_dataset.zip"
    # if output contains an error message, concatenate taxon ID to FINAL_no_records.txt file, otherwise add to successul downloads file, parsedIDs.txt
    #if [[ "\$command" == *"Error: "* ]]; then
    #    echo "$taxonID" >> "FINAL_no_records.txt"
    #else
    #    ...
    #fi


    """

    
}


/*

// OR

process RETRY_FAILED_DOWNLOADS {
    publishDir "${params.outdir}/genome_downloads", mode: params.publish_dir_mode

    input:
    val taxon

    output:
    path "${taxon}_ncbi_dataset", emit: success, optional: true
    path "${taxon}.final_failed", emit: final_failed, optional: true

    script:
    """
    #!/bin/bash
    set -e

    for level in "chromosome" "scaffold" "contig"; do
        if datasets download genome taxon "${taxon}" --assembly-level \$level --dehydrated --filename "${taxon}_ncbi_dataset.zip"; then
            unzip -q "${taxon}_ncbi_dataset.zip"
            datasets rehydrate --directory .
            rm "${taxon}_ncbi_dataset.zip"
            exit 0
        fi
    done

    echo "${taxon}" > "${taxon}.final_failed"
    """
}
*/

