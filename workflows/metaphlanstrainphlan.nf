//

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
include { DOWNLOAD_GENOMES        } from '../subworkflows/local/strain_characterisation/download_genomes.nf'
include { RETRY_DOWNLOADS         } from '../subworkflows/local/strain_characterisation/download_genomes.nf'


include { STRAINPHLAN_PREP_CONSENSUS; STRAINPHLAN_PREP_CLADES  } from '../subworkflows/local/strain_characterisation'
include { STRAINPHLAN_STRAINPHLAN } from '../subworkflows/local/strain_characterisation'




/*
    WORKFLOW: PROFILING & STRAIN_CHARACTERISATION
*/

ch_versions         = Channel.empty()  // Initialise globally
ch_multiqc_files    = Channel.empty() 
ch_final_dbs        = Channel.empty() 

// Define workflow to prep and run MetaPhlAn   
workflow PROFILING {

    take:
    final_input_reads 

    main:

    // Prep specified database or use `metaphlan --install --bowtie2db metaphlan_db_latest`
    if ( params.metaphlan_db ) {
        Channel
            .fromPath( "${params.metaphlan_db}*" )
            .ifEmpty { error "No database files found at ${params.metaphlan_db}*" }
            .toSortedList()
            .set { db_ch }

        UNTAR( db_ch )
        ch_versions = ch_versions.mix(UNTAR.out.versions.first())
        ch_final_dbs = UNTAR.out.untar

    } else if ( params.installdb ) {

        INSTALL_DEPENDENCIES() 
        METAPHLAN_MAKEDB()
        ch_final_dbs = METAPHLAN_MAKEDB.out.db
        ch_versions = ch_versions.mix( METAPHLAN_MAKEDB.out.versions )

        // Define db dir name based on the named index/ db_version if provided, otherwise use metaphlan_db_latest and define this as the directory
        // Set db full path for strainphlan
<<<<<<< HEAD
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
=======
        def db_name = params.metaphlan_index ?: 'metaphlan_db_latest'

        // DEBUG
        println " db_name: ${db_name}      params.metaphlan_index: ${params.metaphlan_index} "
>>>>>>> preprocess_metaphlan

    }

    // Run alignment using MetaPhlAn 
    if ( params.run_metaphlan ) {
        ch_raw_profiles         = Channel.empty()       // Count table/ taxonomy profiles

<<<<<<< HEAD
        // final_input_reads.set{}

        // def (merged_reads, other_reads) = final_input_reads.branch {
        //     merged: it[0].merged
        //     other: true
        // }
        // merged_reads_processed = merged_reads.map { meta, reads -> 
        //     [meta + [single_end: true], reads]
        // }
        // ch_metaphlan_input = merged_reads_processed.mix(other_reads)

=======
>>>>>>> preprocess_metaphlan
        ch_metaphlan_input = final_input_reads.map { meta, reads ->
            if (meta.merged) {
                // For merged paired-end reads, pass only the merged file
                return [meta + [single_end: true, is_merged: true], reads]       // merged file is now technically single-end
            } else {
                // For non-merged reads (both single and paired), pass as is
                return [meta, reads]
            }
        }

        METAPHLAN_METAPHLAN ( 
            ch_metaphlan_input, 
            ch_final_dbs 
        )
        ch_versions        = ch_versions.mix( METAPHLAN_METAPHLAN.out.versions.first() )
        ch_raw_profiles    = ch_raw_profiles.mix( METAPHLAN_METAPHLAN.out.profile )         // Mix profiles for each sample into a single channel

<<<<<<< HEAD
        // Debug - SAM file output: Day_1B: /workdir/e1/4d66d6864353952e072551303fd3e9/Day_1B.sam.bz2
        // METAPHLAN_METAPHLAN.out.sam.view { sample, sam_file ->
        //     println "Debug - initial SAM file output: ${sample.id}: ${sam_file}"
        // } 

        // Ensure SAM file paths are correctly resolved
        // METAPHLAN_METAPHLAN.out.sam
        // .map { id, sam_file ->
        //     // def full_sam_paths = sam_file.collect { sam -> sam.toString() }    // sam_file:workdir, d2, d01f75d86bba50726e630eb671b6cf, 371_BAM134.sam.bz2] // Convert paths to full strings            
        //     // def full_sam_paths = sam_file.collect { it.toString() }      // sam_file:workdir, d2, d01f75d86bba50726e630eb671b6cf, 371_BAM134.sam.bz2] // Convert paths to full strings            
        //     def full_sam_paths = sam_file.startsWith('gs://') ? sam_file : "${params.workdir}/${sam_file}".toString() // TO TRY!!! // Cloud paths remain unaltered or use local paths prefixed with workdir
        //     //def full_sam_paths = sam_file.toString()  //!  [id:371_BAM134, sam_file:/workdir/d2/d01f75d86bba50726e630eb671b6cf/371_BAM134.sam.bz2] // Convert paths to full strings            
        //     [id: id.id, sam_file: full_sam_paths]
        // }
        // .view { println "Debug - Amend SAM file: $it"}
        // .set { ch_sam_files }


=======
        // Debug
        METAPHLAN_METAPHLAN.out.sam.view { sample, sam_file ->
            println "Debug - initial SAM file output: ${sample.id}: ${sam_file}"
        } 
>>>>>>> preprocess_metaphlan

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

    main:

<<<<<<< HEAD
    // strainphlan_db = ch_final_dbs    // mapped and aliased ch_final_dbs to strainphlan_db

    // Convert the database path to a proper channel
    // strainphlan_db = Channel.fromPath( "${params.outdir}/metaphlan_db_*/*.pkl", type: 'file' )
    //     .ifEmpty { "No database found at: ${params.outdir}/metaphlan_db_*/" }
    //     .map { dir -> 
    //         def pkl_path = "${dir}*.pkl"
    //         println "Using database from: ${pkl_path} for srainphlan"
    //         return pkl_path
    //     } // Using database from: /metaphlanstrainphlanTEST/metaphlan_db_latest/*.pkl

=======
>>>>>>> preprocess_metaphlan
    if ( params.run_strainphlan && !params.skip_strainphlan_prep ) {

        if (!params.reference_genomes){
            //downloaded_genomes = Channel.empty() 
            
            //TO BE EDITED BELOW & USE?
            // DOWNLOAD_GENOMES ( clades_list.flatten() )
            // DOWNLOAD_GENOMES.out.downloads_list.view{ "downloaded genomes list ${it}" }

            // download_genomes = DOWNLOAD_GENOMES.out.downloaded_genomes_success
            // reference_genomes = download_genomes

        } else {
<<<<<<< HEAD
            // Use file() to keep original full paths of gc buckets as well as local paths, fromFilePairs() creates tuple id 
            // Channel
            //     .fromFilePairs(
            //         file("${params.reference_genomes}"),     // Match .fna files recursively
            //         checkIfExists: true,
            //         glob: true
            //     )
            //     .map { id, ref_files ->
            //         // def full_paths = ref_files.collect { file -> file.toString() }      // ! Collect paths for all matched files         // [id:GCF_030167985, ref_paths:[/231114-371_517-merged-fastq_ehillman/genome_downloads_complete/33038_ncbi_dataset/data/GCF_030167985.1/GCF_030167985.1_ASM3016798v1_genomic.fna]] 
            //         // def full_paths = ref_files.collect { file -> file.toString().startsWith('gs://') ? file.toString() : "${params.reference_genomes}/**/*.fna" }            // [id:GCF_030167985, ref_paths:[gs://mhra-ngs-dev-ut8t-training/231114-371_517-merged-fastq_ehillman/genome_downloads_complete/33038_ncbi_dataset/**/*.fna]] 
            //         // def full_paths = ref_files.collect { file -> file.startsWith('gs://') ? file.toString() : "${params.reference_genomes}/**/*.fna" }          // [id:GCF_030167985, ref_paths:[gs://mhra-ngs-dev-ut8t-training/231114-371_517-merged-fastq_ehillman/genome_downloads_complete/33038_ncbi_dataset/**/*.fna]] 
            //         // id.original = ref_files.collect { it.toString() } //error original not variabble
            //         def full_paths = ref_files.collect { "${params.reference_genomes}" } 
            //         [id: id, ref_paths: full_paths]
            //     }
            //     .set { all_references }
            // // Debug - Amended full path reference_genomes input: [id:GCA_008121495, ref_paths:[/231114-371_517-merged-fastq_ehillman/genome_downloads_complete/33038_ncbi_dataset/data/GCA_008121495.1/GCA_008121495.1_ASM812149v1_genomic.fna]] 
            // all_references.view { println " Amended full path reference_genomes input: $it " }
=======
            // Keep original full paths of gc buckets as well as local paths
>>>>>>> preprocess_metaphlan
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
                    def strain_id   = row.ref_genome_strainID 
                    def genome_path = row.ref_genome_path ?: error("Missing genome path for strain ${strain_id}")
                    def associated_clade = row.associated_clade
                    def meta = [
                        id: strain_id, 
                        associated_clade: associated_clade
                    ]

                    // Debug
                    // println("Strain ID: ${strain_id}, Genome Path: ${genome_path}")
                    
                    // Return tuple with ID and path
                    [meta, genome_path]        
                }.set { all_references }

            ch_combined_references = all_references
                .map { meta, genome_path -> 
                    genome_path
                }
            ch_combined_references.view { "ch_combined_references: $it" }   
        }


        // Handle missing clades_list otherwise list them and those listed on command line
        if (!clades_list) {
            clades_list = Channel.empty() // Default to an empty channel
        }
        if ( params.strainphlan_clades ) {
            Channel
            .of( "${params.strainphlan_clades}".split(",").toList() )           // .csv file listing clades
            .flatten()
            .set { clade }
        } else {
            clade = clades_list.splitText()
        }

<<<<<<< HEAD
        // clade.view{ "all clades from list: $it" }
=======
        clade.view{ "all clades from list: $it" }
>>>>>>> preprocess_metaphlan

        STRAINPHLAN_PREP_CONSENSUS (
            ch_sam_files,           //     METAPHLAN_METAPHLAN.out.sam,
            ch_final_dbs            //  strainphlan_db, 
        )

        STRAINPHLAN_PREP_CLADES (
            ch_final_dbs,
            clade
        )
        
                
        // Split clade tuple into separate channels
        ch_clade_markers = STRAINPHLAN_PREP_CLADES.out.clade_markers

        ch_clade_markers.map { clade, fna_file ->
            println "Clade: ${clade}, File: ${fna_file}" 
        }

        ch_clade = ch_clade_markers.map { clade, fna_file -> clade }.view { "ch_clade: $it" }        
        ch_marker_file = ch_clade_markers.map { clade, fna_file -> fna_file }.view { "ch_marker_file: $it" } 

        if ( params.strainphlan_db ) {
            strainphlan_db = Channel.fromPath(params.strainphlan_db)
        } else {

            strainphlan_db = Channel.fromPath( "${params.outdir}/metaphlan_db_*", type: 'dir' ) //Channel.fromPath( "${params.outdir}/metaphlan_db_*/**/*.pkl", type: 'file' )
                // .ifEmpty { error "No database found at: ${params.outdir}/metaphlan_db_*/" }
                .map { dir -> 
                    def pkl_path = file(dir).listFiles().findAll { it.name.endsWith('.pkl') }
                    println "Found database .pkl file at: ${pkl_path} for StrainPhlAn"          // Found database .pkl file at: [/metaphlanstrainphlanTEST/metaphlan_db_latest/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl] for StrainPhlAn
                    if (pkl_path.isEmpty()) {
                        error "No .pkl file found in ${dir}"
                    }
                    return pkl_path
                }
                .flatten()
            // strainphlan_db.view {"strainphlan_db: $it "}
        }


        // Combine strain ID with clade for tagging
        // all_references.view { "all_references: $it" }   
  

        STRAINPHLAN_STRAINPHLAN(
            strainphlan_db,                 // First: strainphlan_db
            STRAINPHLAN_PREP_CONSENSUS.out.consensus_markers, // Second: consensus_markers/*.json.bz2
            [],     // ch_combined_references,         // Third: all_references
            ch_profiles,                    // Fourth: merged_profiles
            ch_clade,                       // Fifth: clade
            ch_marker_file                  // Sixth: fna_file
        )
<<<<<<< HEAD
        // // works-sih
=======
        // 
>>>>>>> preprocess_metaphlan
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
        STRAINPHLAN_STRAINPHLAN.out.tre_file.view { ".tre files $it" }

    } else if ( params.skip_strainphlan_prep ) {
        Channel
            .fromPath("${params.outdir}/strainphlan/consensus_markers/*.pkl")
            .ifEmpty { error "No consensus marker files found at ${params.outdir}/strainphlan/consensus_markers/" }
            .set { ch_consensus_markers }

        STRAINPHLAN_STRAINPHLAN (
            ch_consensus_markers,
            ch_final_dbs,
            clade, []
        )

        STRAINPHLAN_STRAINPHLAN.out.tre_file.view { ".tre files $it" }
    }

    // STRAINPHLAN_STRAINPHLAN.out.tre_file.view { ".tre files $it" }


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

    emit:

    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
