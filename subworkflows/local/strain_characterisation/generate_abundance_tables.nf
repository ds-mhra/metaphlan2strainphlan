process GENERATE_ABUNDANCE_TABLES {
    tag "$meta.id"
    publishDir "${params.outdir}/abundance_and_taxons", mode: params.publish_dir_mode

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43ccf4_4' :
        'biocontainers/grep:3.4--hf43ccf4_4' }"

    input:
    tuple val(meta), path(merged_profiles)

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
    echo "Generating species-level abundance table..."
    grep -E "s__" ${merged_profiles} \\
        | grep -v "t__" \\
        | sed "s/^.*|//g" \\
        | cat <(head -n 2 ${merged_profiles}) - \\
        > "merged_abundance_table_species.txt"


    # Generate strain-/ SGB-level abundance table
    echo "Generating strain-/SGB-level abundance table..."
    grep -E "t__" ${merged_profiles} \\
        | sed "s/^.*|//g" \\
        | cat <(head -n 2 ${merged_profiles}) - \\
        | tee "merged_abundance_table_strains.txt" \\
        | grep -E "^t__" | cut -f1 \\
        > "taxons_list.txt" \\
        && echo "Generating taxons list..."


    # Overview includes species name and SGB/clade name along w abundance
    echo "Generating abundance overview..."
    grep -E "t__" ${merged_profiles} \\
        | sed "s/^\\(.*|\\)\\(.*|\\)/\\2/" \\
        | cat <(head -n 2 ${merged_profiles}) - \\
        > "overview_abundance_table_species_strain.txt"

    echo "Abundance tables completed!"

    """
}

