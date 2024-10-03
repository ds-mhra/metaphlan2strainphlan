#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MHRA/metaphlan2strainphlan
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/metaphlanstrainphlan
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { PROFILING               } from './workflows/metaphlanstrainphlan'
// include { STRAIN_CHARACTERISATION } from './workflows/metaphlanstrainphlan'
include { PREPROCESSING                      } from './workflows/preprocessing'
include { PROFILING; STRAIN_CHARACTERISATION } from './workflows/metaphlanstrainphlan'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_metaphlanstrainphlan_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_metaphlanstrainphlan_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_metaphlanstrainphlan_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run metaPhlAn-2-strainPhlAn pipeline
//  
workflow NFCORE_METAPHLANSTRAINPHLAN {
        
    take:
    ch_samplesheet                             // channel: reads from --input_folder or --samplesheet

    main:

    // Define contaminants variable based on params using `file` for file inputs and `null` if none provided
    contaminants = params.shortread_qc_contaminantslist ? file(params.shortread_qc_contaminantslist) : null

    // Run preprocessing
    PREPROCESSING(
        ch_samplesheet,
    )

    // Run profiling with MetaPhlAn
    PROFILING (
        PREPROCESSING.out.final_input_reads,
    )

    // Run characterisation of strains via StrainPhlAn
    STRAIN_CHARACTERISATION (
        PROFILING.out.sam, 
        PROFILING.out.ch_final_dbs,
    )
    
    emit:
    multiqc_report = STRAIN_CHARACTERISATION.out.multiqc_report

} 



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SECOND WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_METAPHLANSTRAINPHLAN (
        PIPELINE_INITIALISATION.out.samplesheet,
//        PIPELINE_INITIALISATION.out.databases,
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_METAPHLANSTRAINPHLAN.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
