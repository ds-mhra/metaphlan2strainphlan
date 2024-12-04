process GENERATE_ABUNDANCE_TABLES {
    tag "$meta.id"
    publishDir "${params.outdir}/abundance_and_taxons", mode: params.publish_dir_mode

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43ccf4_4' :
        'biocontainers/grep:3.4--hf43ccf4_4' }"

    input:
    tuple val(meta), path(merged_profiles)      // metaphlan_results="${params.outdir}/metaphlan/metaphlan_combined_profile_report.txt"

    output:
    tuple val(meta), path("merged_abundance_table_species.txt")     , emit: species_abundance
    tuple val(meta), path("merged_abundance_table_strains.txt")     , emit: strain_abundance
    path "taxons_list.txt"                                          , emit: taxons_list

    script:
    """
    #!/bin/bash
    set -e

    # 1) Generate species-level abundance table
    # 2) Obtain taxon genus and species names only, remove underscore from genus names and add fullstop if any taxon names contains 'sp'
    #   ie isolate headers matching t__ (SGB level), remove excess info and simplify sample names
    # changed to t__ from s_ SO CHECK!!!!!!!!
    echo "Generating species-level abundance table..."
    grep -E "s__" ${merged_profiles} \\
        | grep -v "t__" \\
        | sed "s/^.*|//g" \\
        | tee "merged_abundance_table_species.txt"


    # Generate strain-/ SGB-level abundance table
    echo "Generating strain-/SGB-level abundance table..."
    grep -E "t__" ${merged_profiles} \\
        | sed "s/^.*|//g" \\
        | tee "merged_abundance_table_strains.txt" \\
        | cut -f1 \\
        > "taxons_list.txt" \\
        && echo "Generating taxons list..."


    # Overview includes species name and SGB/clade name along w abundance
    echo "Generating abundance overview..."
    grep -E "t__" ${merged_profiles} \\
        | sed "s/^\\(.*|\\)\\(.*|\\)/\\2/" \\
        > "overview_abundance_table_species_strain.txt"

    echo "Abundance tables completed!"

    """
}


