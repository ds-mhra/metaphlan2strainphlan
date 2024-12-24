/*
    Create series of processes for StrainPhlan
*/

// Runs preparation for StrainPhlan on MetaPhlan's merged profile output
process STRAINPHLAN_PREP_CONSENSUS {

    tag "$meta.id"
    container "quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0"
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"

    input:
    tuple val(meta), path(sam)     // val sam 
    path metaphlan_db_dir          // val ch_final_dbs

    output:
    path "consensus_markers/*", optional: true                  , emit: consensus_markers

    def args = ext.args ?: ''
    """
    #!/bin/bash
    set -e

    PKL_DB=`find -L "${metaphlan_db_dir}" -type f -name "*.pkl"`
    
    # Generate marker files (contains unique marker genes) for each species in the provided database for each sample and store in pickle files (*.pkl)
    printf 'Running sample2markers.py on file: ${sam} using \$PKL_DB'
    mkdir -p consensus_markers

    sample2markers.py -i ${sam} \\
        -o consensus_markers \\
        --database \$PKL_DB \\
        $args \\
        --nproc ${params.cpus}

    """

}

//
// Run ncbi-downloads.sh or other module ?????????????????????
//


process STRAINPHLAN_PREP_CLADES {

    tag "$clade"
    container "quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0"
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"

    input:
    path metaphlan_db_dir
    val clade

    output:
    tuple val(clade), path("clade_markers/*"), optional: true      , emit: clade_markers
    // path "versions.yml"                                                  , emit: versions
    path "rejected_clades.txt", optional: true
    path "parsed_clades.txt"  , optional: true

    def args = ext.args ?: ''
    """
    #!/bin/bash
    set -e

    PKL_DB=`find -L "${metaphlan_db_dir}" -type f -name "*.pkl"`

    # Extract clade markers using *.pkl database from each clade
    mkdir -p clade_markers
    echo "Processing clade: ${clade}"

    # Run extract_markers.py and capture errors
    extract_markers.py --database "\$PKL_DB" --clade "${clade}" --output_dir "clade_markers"

    # If output contains an error message, concatenate to rejected_clades.txt, otherwise successful clades are parsed
    command=\$(extract_markers.py --database "\$PKL_DB" --clade "${clade}" --output_dir "clade_markers" 2>&1)
    if [[ "\$command" == *"Error"* ]]; then
        [[ ! -f "rejected_clades.txt" ]] && touch "rejected_clades.txt"
        echo "${clade}" >> "rejected_clades.txt"
    else
        [[ ! -f "parsed_clades.txt" ]] && touch "parsed_clades.txt" \
        && echo "${clade} processed" 
        echo "${clade}" >> "parsed_clades.txt"
    fi

    """

}


//
// Run StrainPhlan
process STRAINPHLAN_STRAINPHLAN {

    tag "${clade}_with_${strain_id}" 
    container "quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0"
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"

    // Skip process if no data is available for...
    // when:
    // reference_genome && strainphlan_db

    input:
    path strainphlan_db 
    path consensus_markers 
    tuple val(clade), path(reference_genome), path(clade_markers), val(species), val(strain_id) 


            // path strainphlan_db                             // First parameter
            // path consensus_markers                        // Second parameter
            //         // tuple val(meta), path(reference_genome)     //path all_references                // Third parameter
            // // path reference_genome         //all_references
            // tuple val(meta), path(merged_profiles)          // Fourth parameter
            // //tuple val(clade), path(clade_markers)       // Fifth parameter
            // tuple val(clade), path(reference_genome), path(clade_markers)       // Fifth parameter "ADD SPECIES"
            //     // val clade                                 // Fifth parameter
            //     // path clade_fna                            // Sixth parameter
            // path "phylophlan.config"                          // --phylophlan_configuration phylophlan.config 
    


    output:
    tuple val(clade), path("${species}/${strain_id}/*")                                    , emit: strainphlan_outputs
    tuple val(clade), path("${species}/${strain_id}/RAxML_bestTree.*.StrainPhlAn4.tre")    , emit: tre_file        // tuple val(clade), path("${clade}/**/*.tre")     , emit: tre
    // path "versions.yml"                         , emit: versions
    

    script:
    def args = task.ext.args ?: ''

    """
    #!/bin/bash
    set -e

    # Ensure consensus_markers files are present
    echo $consensus_markers

    # Ensure clade_markers directory is not empty, if .fna does not exist for the specified clade, proceed with the next
    if [ -z $clade_markers ]; then
        echo "No clade marker file found for clade, "${clade}". Skipping..."
        exit 0  
    else
        echo "Clade file(s): $clade_markers"
    fi

    # Check reference files exist
    if [ ! -f "$reference_genome" ]; then
        echo "No reference genomes found. Exiting..."
        exit 1
    else
        echo "Reference genomes found: $reference_genome"
    fi

    # Create directories and verify
    mkdir -p "tmp"
    mkdir -p "${species}/${strain_id}"

    echo "Creating directory: ${species}/${strain_id} \n\n"
    if [ ! -d "${species}/${strain_id}" ]; then
        echo "Failed to create directory: ${species}/${strain_id}"
        exit 1
    fi

    # Run StrainPhlAn for each clade and matching reference, while listing the consensus marker files
    strainphlan --samples $consensus_markers \\
        --database $strainphlan_db \\
        --clade_markers $clade_markers \\
        --output_dir "${species}/${strain_id}" \\
        --clade "${clade}" \\
        --phylophlan_mode fast \\
        --mutation_rates --nproc "${params.cpus}" --debug \\
        --tmp "tmp" \\
        $args \\
        --references $reference_genome 

    # Clean up
    echo "Cleaning up temporary files..."
    rm -rf "tmp/*"
    printf '\nStrainPhlAn has finished processing clade: %s, for the reference strain: %s' "${clade}" "${strain_id}"

    """
}


// Add metadata to StrainPhlan tree results
process STRAINPHLAN_METADATA {

    // Container and resources
    container "quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0"
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"

    input:
    tuple val(clade), path(tre)
    path metadata
    // val add_metadata_field

    output:
    tuple val(clade), path("${tre}.*")

    """
    echo "Adding metadata to $tre"

    if [[ -f $metadata ]]; then

        add_metadata_field=\$(head -n 1 $metadata | cut -f 1)
        add_metadata_field == "sampleID"
        
    else
        echo "Error: File not found or file does not have sampleID as the first header."
        exit 1
    fi

    add_metadata_tree.py --ifn_trees ${tre} \\
        --ifn_metadata metadata.txt \\
        --string_to_remove *_merged.fastq.gz \\
        --metadatas ${add_metadata_field}


    plot_tree_graphlan.py \\
        --ifn_tree ${tre}.metadata \\
        --colorized_metadata $colorized_metadata_field \\
        --leaf_marker_size 60 --legend_marker_size 60 \\
        || true     # ignore error, exits with status zero


    """

}

