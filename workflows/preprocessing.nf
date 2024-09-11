/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { UNTAR                  } from '../modules/nf-core/untar/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'                                                                                        
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_metaphlanstrainphlan_pipeline'
include { BBMAP_ALIGN            } from '../modules/nf-core/bbmap/align/main' 
include { BBMAP_BBMERGE          } from '../modules/nf-core/bbmap/bbmerge/main' 
include { BBMAP_BBDUK            } from '../modules/nf-core/bbmap/bbduk/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    contaminants

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()



    // Prep specified database or use `metaphlan --install --bowtie2db metaphlan_db_latest`
    if ( params.db ) {

        process DB_CHECK (
            // Make a channel with the reference database files
            Channel
                .fromPath( "${params.db}*" )
                .ifEmpty { error "No database files found at ${params.db}*" }
                .toSortedList()
                .set { db_ch }

            UNTAR ( db_ch )
            .set { ch_final_dbs }
            ch_versions = ch_versions.mix(UNTAR.out.versions.first())
        )
    } else {
        METAPHLAN_MAKEDB ()
    }

    process INPUT_CHECK (
        
        // Make a channel with reference database files
        Channel
            .fromPath( "${params.db}*" )
            .ifEmpty { error "No database files found at ${params.db}*" }
            .toSortedList()
            .set {db_ch}

        // If the samplesheet exists, convert to tuple/ list by...
        if ( params.samplesheet ){
            // ...creating a channel from the samplesheet, parse it, and branch into single- and paired-end samples
            ch_samplesheet = Channel
                .fromPath(
                    "${params.samplesheet}",
                    checkIfExists: true,
                    glob: true
                )
                .ifEmpty { error "No samplesheet has been added to ${params.samplesheet}*" }   
                .splitCsv(
                    header: true,
                    sep:','
                )
                // Reformat each row of samplesheet so it is easy to pass fastq_1 and fastq_2 columns into tools like fastqc; also handles single-end reads where fastq_2 is null or missing
                .map { row ->
                    tuple(
                        meta = row.sample,
                        [meta, [
                            path(row.fastq_1, checkIfExists: true),
                            row.fastq_2 ? path(row.fastq_2, checkIfExists: true) : null
                        ]]
                    )
                }
                // Branch samples into single- and paired-end sub-channels
                .branch {
                    single: it.fastq_2 == null || it.fastq_2 == ""
                    paired: true
                }
        } else {
        // Make a channel with all of the files from the --input_folder
        inputs = Channel
            .fromFilePairs([
                "${params.input_folder}/*${params.file_spacer}{1,2}${params.file_suffix}"
            ])
            .ifEmpty { error "No file pairs found at ${params.input_folder}/*${params.file_spacer}{1,2}${params.file_suffix}" }
        }
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    

    contaminants = params.shortread_qc_contaminantslist ? file(params.shortread_qc_contaminantslist) : []
    if ( params.shortread_qc_contaminantslist ) {
        if ( params.qc_tool == 'bbduk' && !contaminants.extension.matches(".*(csv|tsv|txt)") ) error "[metaphlanstrainphlanPipeline] ERROR: Contaminants or adapter list requires a different format and extension. Check input: --shortread_qc_contaminantslist ${params.shortread_qc_contaminantslist}"
    }

    // Carry out pre-processing of paired-reads w tools ie bbduk, fastqc 
    if ( !params.skip_preprocessing_qc ) {
        if ( params.perform_shortread_qc && params.qc_tool == 'bbduk' ) {
            ch_shortreads_preprocessed = BBMAP_BBDUK ( 
                ch_samplesheet.paired, contaminants
                )
            ch_versions = ch_versions.mix( BBMAP_BBDUK.out.versions )
        } else {
            ch_shortreads_preprocessed = ch_samplesheet.paired # INPUT_CHECK.out.paired
        }

        //
        // MODULE: Run FastQC
        //
        if ( params.perform_shortread_qc && params.qc_tool == 'fastqc' ) {    
            FASTQC ( 
                ch_shortreads_preprocessed  

            )
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        }

        //
        // MODULE: Run BBmerge
        //
        if (params.bbmerge_pairs ) {
            ch_merged_reads = BBMAP_BBMERGE (
                ch_shortreads_preprocessed
            )
            ch_versions = ch_versions.mix(BBMAP_BBMERGE.out.versions.first())
            // If not merging runs, use bbmerged reads as final input for metaphlan
            if ( !params.perform_runmerging ) {
                final_input_reads = ch_merged_reads
            }
        }
    }
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

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
