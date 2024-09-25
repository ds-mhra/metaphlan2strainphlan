#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

/*
    Create series of processes for StrainPhlan
*/


// Runs preparation for StrainPhlan on MetaPhlans merged profile output
process STRAINPHLAN_PREP {

    // Container to run process
    container "quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0"

    // Resources used
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"
    publishDir "${params.outdir}/strainphlan", "${params.publish_dir_mode}", overwrite: true

    input:
    path sam 
    path ch_final_dbs
    val clade

    output:
    path "consensus_markers/*", optional: true
    path "clade_markers/*", optional: true


    """#!/bin/bash
    set -e
    
    # Create required directories
    mkdir -p {consensus_markers,clade_markers}

    # Generate marker files (contains unique marker genes) for each species in the provided database for each sample and store in pickle files (*.pkl)
    printf 'Running sample2markers.py on file: %s\n' "${sam}"
    
    sample2markers.py -i ${sam} \
        -o consensus_markers \
        --database "${params.strainphlan_db}" \
        --nproc "${params.nproc}" \
        $args

    # Extract clade markers using *.pkl database from each clade
    extract_markers.py --database "${params.strainphlan_db}" \
        --clade ${clade} \
        --output_dir "clade_markers"

    """

}


//
// Run StrainPhlan
process STRAINPHLAN_STRAINPHLAN {

    // Container to run process
    container "quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0"

    // Resources used
    cpus "${params.cpus}"
    memory "${params.memory_gb}.GB"
    publishDir "${params.outdir}/strainphlan/strainphlan_output", "${params.publish_dir_mode}", overwrite: true

    input:
    path ch_final_dbs
    path "consensus_markers/*.pkl"
    tuple val(clade), path("clade_markers/*.fna")
    //val clade
    //path "clade_markers/*"
    path "phylophlan.config"

    output:
    path "${clade}/*", emit: all
    tuple val(clade), path("strainphlan_output/${clade}/*.tre"), emit: tre


    """#!/bin/bash
    set -e
    
    mkdir -p ${params.outdir}/{tmp,strainphlan_output/${clade}}

    printf '\n\nClade: %s\n' "$clade"
    strainphlan --samples consensus_markers/*.pkl \
            --database "${params.strainphlan_db}"/*.pkl \
            --clade_markers "${params.outdir}/clade_markers/${clade}.fna" \
            --output_dir "${params.outdir}/strainphlan_output/${clade}" \
            --clade "${clade}" \
            --phylophlan_mode fast \
            --phylophlan_configuration phylophlan.config \
            --mutation_rates --nproc "${params.nproc}" --debug \
            --tmp "${params.outdir}/tmp" \
            $args
    printf '\n\nStrainPhlAn has finished processing clade: %s\n' "$clade"

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
    tuple val(clade), path(tre), path("metadata.txt")

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
            add_metadata_tree.py --ifn_trees $tre_path \
                --ifn_metadata metadata.txt \
                --string_to_remove *_merged.fastq.gz \
                --metadatas sampleID
        done
    done

    """
}

