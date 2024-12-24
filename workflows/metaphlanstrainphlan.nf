/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_metaphlanstrainphlan_pipeline'
include { INSTALL_DEPENDENCIES   } from "../modules/local/installing_dependencies.nf"
include { BLAST_MAKEBLASTDB      } from '../modules/nf-core/blast/makeblastdb/main'
include { BOWTIE2_ALIGN          } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD          } from '../modules/nf-core/bowtie2/build/main'
include { METAPHLAN_MAKEDB       } from '../modules/nf-core/metaphlan/makedb/main'                                        
include { METAPHLAN_METAPHLAN    } from '../modules/nf-core/metaphlan/metaphlan/main'                                                            
include { METAPHLAN_MERGEMETAPHLANTABLES } from '../modules/nf-core/metaphlan/mergemetaphlantables/main'

include { GENERATE_ABUNDANCE_TABLES } from '../subworkflows/local/strain_characterisation/generate_abundance_tables.nf'
// include { DOWNLOAD_GENOMES        } from '../subworkflows/local/strain_characterisation/download_genomes.nf'
// include { RETRY_DOWNLOADS         } from '../subworkflows/local/strain_characterisation/download_genomes.nf'

include { STRAINPHLAN_PREP_CONSENSUS; STRAINPHLAN_PREP_CLADES  } from '../subworkflows/local/strain_characterisation'
include { STRAINPHLAN_STRAINPHLAN; STRAINPHLAN_METADATA        } from '../subworkflows/local/strain_characterisation'




/*
    WORKFLOW: PROFILING & STRAIN_CHARACTERISATION
*/

ch_final_dbs        = Channel.empty() 

// Define workflow to prep and run MetaPhlAn   
workflow PROFILING {

    take:
    final_input_reads
    ch_versions

    main:

    // Prep specified database or use `metaphlan --install --bowtie2db metaphlan_db_latest`
    if ( params.metaphlan_db ) {
        Channel
            .fromPath( "${params.metaphlan_db}" )
            .ifEmpty { error "No database found at ${params.metaphlan_db}." }
            // .toSortedList() //.toList()
            .set { ch_final_dbs }

        // UNTAR( db_ch )
        // ch_versions = ch_versions.mix(UNTAR.out.versions.first())
        // ch_final_dbs = UNTAR.out.untar

    } else if ( params.installdb ) {

        INSTALL_DEPENDENCIES() 
        METAPHLAN_MAKEDB()
        ch_final_dbs = METAPHLAN_MAKEDB.out.db
        ch_versions = ch_versions.mix( METAPHLAN_MAKEDB.out.versions )

        // Define db dir name based on the named index/ db_version if provided, otherwise use metaphlan_db_latest and define this as the directory
        // Set db full path for strainphlan
        // def db_name = params.metaphlan_index ? "metaphlan_db_${params.metaphlan_index}" : 'metaphlan_db_latest'
        // OR simply: 
        def db_name = params.metaphlan_index ?: 'metaphlan_db_latest'

        // def strainphlan_db = "${db_name}/*.pkl"          // must add tp println      strainphlan_db: ${strainphlan_db}
        // def strainphlan_db = db_name.findAll { it.endsWith('.pkl') }


        //params.metaphlan_db = db_name  
        //params.strainphlan_db = "${params.outdir}/${params.metaphlan_index}/*.pkl"

        // DEBUG
        println " db_name: ${db_name}      params.metaphlan_index: ${params.metaphlan_index} "
        //      println " params.metaphlan_index: ${params.metaphlan_index}    db_name: ${db_name}      params.metaphlan_index: ${params.metaphlan_index}      params.strainphlan_db: ${params.strainphlan_db}"
        //       params.metaphlan_db: null   params.metaphlan_index: mpa_vOct22_CHOCOPhlAnSGB_202212         params.strainphlan_db: gs://mhra-ngs-dev-ut8t-training/metaphlanstrainphlanTEST/mpa_vOct22_CHOCOPhlAnSGB_202212/*.pkl

    }

    // Run alignment using MetaPhlAn 
    if ( params.run_metaphlan ) {
        ch_raw_profiles         = Channel.empty()       // Count table/ taxonomy profiles

        // final_input_reads.set{}

        // def (merged_reads, other_reads) = final_input_reads.branch {
        //     merged: it[0].merged
        //     other: true
        // }
        // merged_reads_processed = merged_reads.map { meta, reads -> 
        //     [meta + [single_end: true], reads]
        // }
        // ch_metaphlan_input = merged_reads_processed.mix(other_reads)

        ch_metaphlan_input = final_input_reads.map { meta, reads ->
            if (meta.merged) {
                // For merged paired-end reads, pass only the merged file
                return [meta + [single_end: true, is_merged: true], reads]       // merged file is now technically single-end
            } else {
                // For non-merged reads (both single and paired), pass as is
                return [meta, reads]
            }
        }
        ch_metaphlan_input.view{ "ch_metaphlan_input $it " }

        METAPHLAN_METAPHLAN ( 
            ch_metaphlan_input, 
            ch_final_dbs 
        )
        ch_versions        = ch_versions.mix( METAPHLAN_METAPHLAN.out.versions.first() )
        ch_raw_profiles    = ch_raw_profiles.mix( METAPHLAN_METAPHLAN.out.profile )         // Mix profiles for each sample into a single channel

        // Debug
        // METAPHLAN_METAPHLAN.out.sam.view { sample, sam_file ->
        //     println "Debug - initial SAM file output: ${sample.id}: ${sam_file}"
        // }


        // Merge all MetaPhlAn profiles
        // First re-map each sample id and its profile, then group all profiles by sample name and merge
        METAPHLAN_MERGEMETAPHLANTABLES ( 
            ch_raw_profiles.map{ [ [id:'all_samples'], it[1] ] }.groupTuple(sort: { it.getName() })
        )
        
        ch_versions        = ch_versions.mix( METAPHLAN_MERGEMETAPHLANTABLES.out.versions.first() )
    }

    emit:
    ch_sam_files = METAPHLAN_METAPHLAN.out.sam
    ch_raw_profiles
    ch_profiles = METAPHLAN_MERGEMETAPHLANTABLES.out.txt       // ERROR ~ Cannot emit a multi-channel output: ch_profiles
    ch_final_dbs
    versions = ch_versions

}


// Extracts clades from metaphlan results for input into strainphlan
workflow METAPHLAN_TO_STRAINPHLAN {

    take:
    ch_profiles

    main:
    GENERATE_ABUNDANCE_TABLES (ch_profiles)

    emit:
    clades_list         = GENERATE_ABUNDANCE_TABLES.out.taxons_list

}


// Define workflow to prep and/or run StrainPhlAn for strain characterisation
workflow STRAIN_CHARACTERISATION {

    take:
    ch_sam_files 
    ch_final_dbs    // ${params.strainphlan_db}   
    clades_list     // TBD
    ch_profiles     // Added
    ch_multiqc_files
    ch_versions


    main:

    // strainphlan_db = ch_final_dbs    // mapped and aliased ch_final_dbs to strainphlan_db

    // Convert the database path to a proper channel
    // strainphlan_db = Channel.fromPath( "${params.outdir}/metaphlan_db_*/*.pkl", type: 'file' )
    //     .ifEmpty { "No database found at: ${params.outdir}/metaphlan_db_*/" }
    //     .map { dir -> 
    //         def pkl_path = "${dir}*.pkl"
    //         println "Using database from: ${pkl_path} for srainphlan"
    //         return pkl_path
    //     } // Using database from: /metaphlanstrainphlanTEST/metaphlan_db_latest/*.pkl

    if (!params.reference_genomes){
        //downloaded_genomes = Channel.empty() 
        
        //TO BE EDITED BELOW & USE?
        // DOWNLOAD_GENOMES ( clades_list.flatten() )
        // DOWNLOAD_GENOMES.out.downloads_list.view{ "downloaded genomes list ${it}" }

        // download_genomes = DOWNLOAD_GENOMES.out.downloaded_genomes_success
        // reference_genomes = download_genomes

    } else {
        Channel
            .fromPath(
                params.reference_genomes,
                checkIfExists: true,
                glob: true
            )
            .ifEmpty { error "No genome list has been found at ${params.reference_genomes}" } 
            .splitCsv(
                header: true,
                sep:','
            )
            // Reformat each row of samplesheet so it is easy to pass fastq_1 and fastq_2 columns into tools like fastqc; also handles single-end reads where fastq_2 is null or missing
            .map { row ->

                // Map sample data with file checks by creating a tuple of sample metadata (id) and fastq_files
                def species_name     = row.species_name
                def strain_id        = row.ref_genome_strainID 
                def genome_path      = row.ref_genome_path ?: error("Missing genome path for strain ${strain_id}")
                def associated_clade = row.associated_clade ?: error("Missing associated clade for the species ${species_name}. Add clade name to ${params.reference_genomes}.")
                def meta = [
                    species_name: species_name,
                    strain_id: strain_id, 
                    associated_clade: associated_clade
                ]

                // Debug
                // println("Species: ${species_name}, Strain ID: ${strain_id}, Genome Path: ${genome_path}")
                
                // Return tuple with ID and path
                [meta, genome_path]        
            }.set { all_references }

        // all_references.view { "all_references: $it" }

    }

    if ( params.run_strainphlan && !params.skip_strainphlan_prep ) {

        // Create channel for clades listed on command line, otherwise use clades_list generated from abundance table
        if ( params.strainphlan_clades ) {
            Channel
            .of( "${params.strainphlan_clades}".split(",").toList() )           // comma separated list of clades
            .flatten()
            .set { clade }
        } else if (clades_list) {
            clade = clades_list.splitText().toList()
        } else {
            clade = Channel.empty() // Default to an empty channel
        }
        // AKA
        // clade = params.strainphlan_clades 
        //     ? Channel.of("${params.strainphlan_clades}".split(",").toList().flatten() )
        //     : clades_list 
        //         ? clades_list.splitText().toList()
        //         : Channel.empty()       

        // clade.view{ "all clades from list: $it" }

        STRAINPHLAN_PREP_CONSENSUS (
            ch_sam_files,           //     METAPHLAN_METAPHLAN.out.sam,
            ch_final_dbs            //  strainphlan_db, 
        )

        STRAINPHLAN_PREP_CLADES (
            ch_final_dbs,
            clade
        )
        
                
        // Split clade tuple into separate channels
        STRAINPHLAN_PREP_CLADES.out.clade_markers
            .set { ch_clade_markers }
        ch_clade_markers.view { " ch_clade_markers: $it "}  // [t__SGB4936, /workdir/68/6905dbe8e85007019549cf74e17572/clade_markers/t__SGB4936.fna]

            // ch_clade_markers.map { clade, fna_file ->
            //     println "Clade: ${clade}, File: ${fna_file}" 
            // }
            //ch_clade = ch_clade_markers.map { clade, fna_file -> clade }    //.view { "ch_clade: $it" }        
            //ch_marker_file = ch_clade_markers.map { clade, fna_file -> fna_file }   //.view { "ch_marker_file: $it" } 
            
            // ch_combined_references = ch_clade_markers
            //     .map { clade, fna_file ->
            //         // Read genome_list.csv and find matching references for the clade
            //         def matched_rows = params.reference_genomes.findAll { row -> row.associated_clade == clade }
            //         if (matched_rows.isEmpty()) {
            //             error "No matching reference genome found for clade ${clade} in ${params.reference_genomes}"
            //         }
            //         // Collect reference genome paths
            //         def genome_paths = matched_rows.collect { row -> row.ref_genome_path }
            //         [clade, genome_paths, fna_file]
            //     }
            //     .view { "Clade: $it[0], Reference Genomes: $it[1], Marker File: $it[2]" }

        
        // all_references.view { " all_references: $it " }   // [[species_name:Roseburia hominis, strain_id:GCF_902387955.1, associated_clade:t__SGB4936], gs://bucket_name/././Roseburia hominis_ncbi_dataset/data/GCF_902387955.1/GCF_902387955.1_UHGG_MGYG-HGUT-02517_genomic.fna]
        
        
        // Modify the clade matching process
        ch_combined_references = ch_clade_markers
            .combine(all_references)
            .filter { clade, fna_marker_file, meta, genome_path -> 
                boolean matched = meta.associated_clade == clade 
                if (!matched) {
                    println "No match for Clade ${clade}. Meta clade: ${meta.associated_clade}"
                }
                matched
            }
            .map { clade, fna_marker_file, meta, genome_path ->
                def species_clean = meta.species_name.replace(' ', '_')         // Replace special characters in species name
                [
                // tuple(      
                    clade,
                    genome_path,
                    fna_marker_file,
                    species_clean,
                    meta.strain_id 
                // )
                ]
            }
            // ch_combined_references.view { tuple -> 
            //     "Clade-Reference Combination: Tuple: ${tuple}" 
            // } 
            // Clade-Reference Combining: Tuple: [t__SGB4936, gs://mhra-ngs-dev-ut8t-training/231114-371_517-merged-fastq_ehillman/genome_downloads_complete/Roseburia hominis_ncbi_dataset/data/GCF_000225345.1/GCF_000225345.1_ASM22534v1_genomic.fna, /workdir/68/6905dbe8e85007019549cf74e17572/clade_markers/t__SGB4936.fna, Roseburia_hominis, GCF_000225345.1]
            // Clade-Reference Combining: Tuple: [t__SGB4936, gs://mhra-ngs-dev-ut8t-training/231114-371_517-merged-fastq_ehillman/genome_downloads_complete/Roseburia hominis_ncbi_dataset/data/GCF_028743455.1/GCF_028743455.1_ASM2874345v1_genomic.fna, /workdir/68/6905dbe8e85007019549cf74e17572/clade_markers/t__SGB4936.fna, Roseburia_hominis, GCF_028743455.1]
            // Clade-Reference Combining: Tuple: [t__SGB4936, gs://mhra-ngs-dev-ut8t-training/231114-371_517-merged-fastq_ehillman/genome_downloads_complete/Roseburia hominis_ncbi_dataset/data/GCF_902387955.1/GCF_902387955.1_UHGG_MGYG-HGUT-02517_genomic.fna, /workdir/68/6905dbe8e85007019549cf74e17572/clade_markers/t__SGB4936.fna, Roseburia_hominis, GCF_902387955.1]

            // OR
            // ch_combined_references.view { clade, genome, marker_file, species, strain_id -> 
            //     "Clade-Reference Matching: \n  Clade: ${clade}   \n  Species: ${species}    \n Genome: ${genome} \n  Marker: ${marker_file}" 
            // }



            // ch_combined_references2 = ch_combined_references
            //     .map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) }
            //     .toList()
            //     .view { "ch_combined_references2222 $it" } // Ensures each tuple is separate



                                // ch_clade_markers
                                    // .map { clade, fna_file ->
                                    // def clade_compare = clade.toString()
                                    // all_references
                                    //     .filter { meta, genome_path -> 
                                    //         println "  2 Comparing: ${meta.associated_clade} vs ${clade_compare}" //Comparing: t__SGB4936 vs DataflowBroadcast around DataflowStream[?]
                                    //         meta.associated_clade == clade_compare
                                    //     }
                                    //     .map { meta, genome_path -> 
                                            // [clade_compare, genome_path, file(fna_file)] // Create the tuple
                                    //     }
                                    // }
                                    // .set { ch_combined_references }

                                        // ch_combined_references.view { "ch_combined_references $it" }

                                // Modify the clade matching process
                                // ch_combined_references = ch_clade_markers
                                //     .map { clade, fna_file ->
                                //         // Read reference genomes and find matching references for the clade
                                //         def matched_references = all_references
                                //             .filter { meta, genome_path -> 
                                //                 println "Comparing: ${meta.associated_clade} vs ${clade}"
                                //                 meta.associated_clade == clade }     // filter meta and genome paths, ensure 'associated_clade' has the variable 'clade'
                                //             .map { meta, genome_path -> genome_path }
                                //             .toList()
                                //             .get()

                                //         if (matched_references.isEmpty()) {
                                //             error "No matching reference genome found for clade ${clade} in ${params.reference_genomes}"
                                //         }

                                //         [clade, matched_references, file(fna_file)]
                                //     }.view { "Clade: $it[0], Reference Genomes: $it[1], Marker File: $it[2]" }



        if ( params.strainphlan_db ) {
            strainphlan_db = Channel.fromPath(params.strainphlan_db)
        } else {

            strainphlan_db = ch_final_dbs
                    // .ifEmpty { error "No database found at: ${params.outdir}/metaphlan_db_*/" }          // causes error at the start of pipeline
                .map { dir -> 
                    def pkl_path = file(dir).listFiles().findAll { it.name.endsWith('.pkl') }
                    println "Found database .pkl file at: ${pkl_path} for StrainPhlAn"          // 
                    if (pkl_path.isEmpty()) {
                        error "No .pkl file found in ${dir}"
                    }
                    return pkl_path
                }
                .flatten() 
        }
        // strainphlan_db.view {"strainphlan_db: $it "}

        STRAINPHLAN_STRAINPHLAN(
            strainphlan_db.first(),  // since channel contains 1 entry, .first() converts channel into a value channel, allowing it to be reused, identify entry w `strainphlan_db.count().view { "No of database entries: $it" } `
            STRAINPHLAN_PREP_CONSENSUS.out.consensus_markers.collect(),
            ch_combined_references 
        )
        

        // // works-ish
        // STRAINPHLAN_STRAINPHLAN (
        //     STRAINPHLAN_PREP_CONSENSUS.out.consensus_markers,
        //     strainphlan_db, //ch_final_dbs, 
        //     all_references,
        //     ch_profiles,
        //     ch_clade,
        //     ch_marker_file
            
        //     // []
        //     // STRAINPHLAN_PREP_CLADES.out.clade_markers
        //     // clade
        // )
        // STRAINPHLAN_STRAINPHLAN.out.tre_file.view { ".tre files $it" }


        if ( params.metadata ) {
            ch_metadata = Channel.fromPath("${params.metadata}", checkIfExists: true)
                .ifEmpty { error "No metadata file has been found at: ${params.metadata}." }

            STRAINPHLAN_METADATA(
                STRAINPHLAN_STRAINPHLAN.out.tre_file,
                ch_metadata
            )
        }

    } else if ( params.skip_strainphlan_prep ) {
        Channel
            .fromPath("${params.consensus_markers}/*")   // !! OR "${params.outdir}/**/*.pkl"
            .ifEmpty { error "No consensus marker files (.pkl or .json.bz2) found at ${params.outdir}/strainphlan/consensus_markers/" }
            .set { ch_consensus_markers }

        Channel
            .fromPath("${params.clade_markers}/*.fna")
            .ifEmpty { error "No clade marker files (.fna) found at ${params.clade_markers}." }
            .set { ch_clade_markers }
        
        // Create clade channel
        clade = params.strainphlan_clades 
            ? Channel.of("${params.strainphlan_clades}".split(",").toList().flatten() )
            : clades_list 
                ? clades_list.splitText().toList()
                : Channel.empty() 

        // Match clade channel to the clade in the reference_genome file and remap into tuple
        ch_combined_references = ch_clade_markers
            .combine(all_references)
            .filter { clade, fna_marker_file, meta, genome_path -> 
                boolean matched = meta.associated_clade == clade 
                if (!matched) {
                    println "No match for StrainPhlAn clade input: ${clade}  vs  ${params.reference_genomes}'s associated clade(s) ${meta.associated_clade}"
                }
                matched
            }
            .map { clade, fna_marker_file, meta, genome_path ->
                def species_clean = meta.species_name.replace(' ', '_')         // Replace special characters in species name
                [
                    clade,
                    genome_path,
                    fna_marker_file,
                    species_clean,
                    meta.strain_id 
                ]
            }

        STRAINPHLAN_STRAINPHLAN (
            strainphlan_db.first(),       
            ch_consensus_markers.collect(),
            ch_combined_references
        )
        // STRAINPHLAN_STRAINPHLAN.out.tre_file.view { ".tre files $it" }

        STRAINPHLAN_METADATA(
            STRAINPHLAN_STRAINPHLAN.out.tre_file,
            ch_metadata
        )
    }


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
    

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    // DEBUG
    // ch_multiqc_files.subscribe { println "Files passed to MultiQC: $it" }


    emit:

    multiqc_report = MULTIQC.out.report.toList() 
    versions       = ch_versions                 
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
