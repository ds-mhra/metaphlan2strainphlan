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
    path "consensus_markers/*", optional: true                     , emit: consensus_markers

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
    //publishDir "${params.outdir}/strainphlan", "${params.publish_dir_mode}", overwrite: true

    input:
    path metaphlan_db_dir          // val ch_final_dbs
    val clade

    output:
    tuple val(clade), path("clade_markers/${clade}.fna")  , optional: true      , emit: clade_markers     // changed from path "clade_markers/*" to tuple val(clade), path("clade_markers/*")

    def args = ext.args ?: ''
    """
    #!/bin/bash
    set -e

    PKL_DB=`find -L "${metaphlan_db_dir}" -type f -name "*.pkl"`

    # Extract clade markers using *.pkl database from each clade
    mkdir -p clade_markers
    echo "Processing clade: ${clade}"

    # Run extract_markers.py and capture errors
    if ! extract_markers.py --database "\$PKL_DB" --clade "${clade}" --output_dir "clade_markers"; then
        echo "Error: Clade ${clade} failed. Logging to rejected_clades.txt."
        echo "${clade}" >> rejected_clades.txt
        exit 0  # Allow the process to continue with other clades
    fi

    """

}


//
// Run StrainPhlan
process STRAINPHLAN_STRAINPHLAN {

    tag "$meta"
    // tag "${meta.id}_${meta.clade}"
    container "quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0"
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"
    // publishDir "${params.outdir}/strainphlan/", "${params.publish_dir_mode}", overwrite: true

    input:
    path strainphlan_db                // First parameter
    path "consensus_markers/*"         // Second parameter
    // tuple val(meta), path(reference_genome)     //path all_references                // Third parameter
    path reference_genome
    tuple val(meta), path(merged_profiles) // Fourth parameter
    val clade                          // Fifth parameter
    val fna_file                       // Sixth parameter
    //path "phylophlan.config"                          // --phylophlan_configuration phylophlan.config 
    

    output:
    // val "${clade}/*"                                               , emit: all
    // path "*.tre", emit: tre_file
    // path "${clade}/*"                           , emit: all
    tuple val(clade), path("${clade}/*.tre")    , emit: tre_file        // tuple val(clade), path("${clade}/**/*.tre")     , emit: tre
    path "species_to_clade_mapping.txt"

    script:
    def args = task.ext.args ?: ''

    """
    #!/bin/bash
    set -e

    # Ensure consensus_markers directory is not empty
    if [ ! "\$(ls -A "consensus_markers/*" 2>/dev/null)" ]; then
        echo "No consensus markers found. Exiting..."
        exit 1
    fi

    # Ensure clade_markers directory is not empty or .fna does not exist for the specified clade, proceed with the next
    if [ ! "\$(ls -A "${fna_file}" 2>/dev/null)" ]; then
        echo "No clade marker file found for clade, "${clade}". Skipping..."
        exit 0  
    fi

    # Check reference files exist
    if [ ! "\$(ls -A ${reference_genome} 2>/dev/null)" ]; then
        echo "No reference genomes found. Exiting..."
        exit 1
    fi


    # Prep and create mapping file
    mkdir -p "{tmp,\${clade}}"
    mapping_file="species_to_clade_mapping.txt"
    echo "species_name clade" > "\$mapping_file"
    grep -E "t__" ${merged_profiles} \\
        | awk -F'|' '{print \$1,\$2}' \\
        | sed 's/s__//g' \\
        >> "\$mapping_file"


    # Run StrainPhlAn loop for each clade and matching reference
    for fna in ${reference_genome}
    do
        echo "Processing reference genome: \$fna"
                #ref_name=\$(basename \$fna .fna)
                #output_dir="strainphlan_output/\${clade}/\${ref_name}"


        # Attempt multiple name matching
        ref_match=\$(basename "\${fna}" .fna)
        clade_match=\$(grep "\$ref_match" "\$mapping_file" | cut -d' ' -f2)

        if [ -z "\$clade_match" ]; then
            echo "No matching clade for \${fna}. Skipping..."
            continue
        fi


        printf '\n\nProcessing clade: %s\n ' "\$clade_match"
        printf '\n\nProcessing reference genome: %s\n ' "\$ref_match"

        # StrainPhlAn 
        strainphlan --samples "consensus_markers/*.json.bz2" \\
            --database "${strainphlan_db}" \\
            --clade_markers "${fna_file}" \\
            --output_dir "${clade}" \\
            --clade "${clade}" \\
            --phylophlan_mode fast \\
            --mutation_rates --nproc "${params.cpus}" --debug \\
            --tmp "tmp" \\
            $args \\
            --references "\${fna}" 

    done

    # Clean up
    echo "Cleaning up temporary files..."
    rm -rf tmp/*
    printf '\n\nStrainPhlAn has finished processing clade: %s\n' "\$clade"

    """
}


// Add metadata to StrainPhlan tree results
process STRAINPHLAN_METADATA {

    // Container to run process
    container "quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0"

    // Resources used
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"
    publishDir "${params.outdir}/strainphlan/strainphlan_output", "${params.publish_dir_mode}", overwrite: true, pattern: "*.tre.metadata"

    input:
    tuple val(clade), path(tre)
    path("metadata.txt")

    output:
    tuple val(clade), path("${tre}.*")

    """
    echo "Adding metadata to ${tre}"

    for clade_path in "${params.outdir}"/strainphlan/strainphlan_output/*
    do
        clade=\$(basename "$clade_path")
        printf '%s\n' "Strain/ Clade: $clade"
        ls "${clade_path}/RAxML_bestTree.${clade}.StrainPhlAn4.tre" | while IFS= read -r tre_path
        do
            treeFile=\$(basename "$tre_path")
            printf '%s\n' "$treeFile"
            add_metadata_tree.py --ifn_trees $tre_path \\
                --ifn_metadata metadata.txt \\
                --string_to_remove *_merged.fastq.gz \\
                --metadatas sampleID
        done
    done

    """
}

