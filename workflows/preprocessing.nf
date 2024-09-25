/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UNTAR                  } from '../modules/nf-core/untar/main'
include { INPUT_CHECK            } from '../subworkflows/input_check'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_PRE   } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_POST  } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_metaphlanstrainphlan_pipeline'
include { BBMAP_ALIGN            } from '../modules/nf-core/bbmap/align/main' 
include { BBMAP_BBDUK as BBMAP_BBDUK_PAIRED; BBMAP_BBDUK as BBMAP_BBDUK_SINGLE } from '../modules/nf-core/bbmap/bbduk/main'

// Odd issue with error: `Process 'BBMAP_BBMERGE' has been already used -- If you need to reuse the same component, include it with a different name or include it in a different workflow context`
include { BBMAP_BBMERGE as BBMAP_BBMERGE_PRE   } from '../modules/nf-core/bbmap/bbmerge/main' 


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
    WORKFLOW: PREPROCESSING
*/

workflow PREPROCESSING {

    take:
    ch_samplesheet     // channel: samplesheet read in from --input 
    
    main:
    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()

    // If the samplesheet exists, convert to tuple/ list by...
    if ( params.samplesheet ){

/*
        // ...verify samplesheet and...
        INPUT_CHECK(
            params.samplesheet
        )
        .set { ch_verified_samplesheet }
        .view()      // checks by viewing samplesheet
*/

        // ...create a channel from the samplesheet, parse it as separated R1nR2 reads, and branch into single- and paired-end samples
        ch_verified_samplesheet = Channel
            .fromPath(
                "${params.samplesheet}",
                checkIfExists: true,
                glob: true
            )
            .ifEmpty { error "No samplesheet has been found at ${params.samplesheet}" } 
            .splitCsv(
                header: true,
                sep:','
            )
            // Reformat each row of samplesheet so it is easy to pass fastq_1 and fastq_2 columns into tools like fastqc; also handles single-end reads where fastq_2 is null or missing
            .map { row ->
                tuple(
                    row.sample,     // sample metadata
                    [
                        file(row.fastq_1, checkIfExists: true),
                        row.fastq_2 ? path(row.fastq_2, checkIfExists: true) : null
                    ]
                )
            }
            // Branch samples into single- and paired-ends sub-channels from mapped tuple(sample, [fastq_1, fastq_2]) - second object within the second element of the array (i.e. fastq_2) is null == single-ends 
            .branch {
                single: it[1][1] == null
                paired: it[1][1] != null
            }
        // Assign the resulting channels from the branch
        ch_samplesheet_single = ch_verified_samplesheet.single
        ch_samplesheet_paired = ch_verified_samplesheet.paired

    } else {
    // Make a channel with all of the files from the --input_folder
    ch_input_files = Channel
        .fromFilePairs([
            "${params.input_folder}/*${params.file_spacer}{1,2}${params.file_suffix}"
        ])
        .ifEmpty { error "No file pairs found at ${params.input_folder}/*${params.file_spacer}{1,2}${params.file_suffix}" }
    }


    //
    // SUB-WORKFLOW: QC checks for adapters, contaminants etc as well as preprocessing of single- paired-end reads w tools ie bbduk, fastqc, fastp

    //
    // MODULE: Run FastQC before processing
    //
    FASTQC_PRE ( 
        ch_samplesheet_paired.mix(ch_samplesheet_single) 
    )
    ch_versions = ch_versions.mix(FASTQC_PRE.out.versions.first())
    

    // Perform QC if either perform_shortread_qc or a qc_tool is selected
    if (params.perform_shortread_qc || params.qc_tool) {
        
        // Split tools and iterate over them
        params.qc_tool.split(',').each { tool ->

            // Trim the tool to remove any spaces
            tool = tool.trim()

            // Error message
            if (!['bbduk', 'fastp', 'fastqc', 'bbmerge'].contains(tool)) {
                error "[metaphlan2strainphlanPipeline] ERROR: Unsupported qc_tool: ${tool}"
            }

            // If the qc_tool is 'bbduk', run BBDUK for QC
            if (tool == 'bbduk') {

                contaminants = params.shortread_qc_contaminantslist ? file(params.shortread_qc_contaminantslist) : []

                if (params.shortread_qc_contaminantslist && !contaminants.extension.matches(".*(csv|tsv|txt)")) {
                    error "[metaphlan2strainphlanPipeline] ERROR: Contaminants or adapter list requires a different format and/or extension. Check input: --shortread_qc_contaminantslist ${params.shortread_qc_contaminantslist}"
                }

                //
                // MODULE: Run BBMAP_BBDUK
                // Trim and filter paired-end reads with BBDUK
                ch_shortreads_paired_preprocessed = BBMAP_BBDUK_PAIRED ( 
                    ch_samplesheet_paired, contaminants 
                ).reads

                // Trim and filter single-end reads with BBDUK
                ch_shortreads_single_preprocessed = BBMAP_BBDUK_SINGLE (
                    ch_samplesheet_single, contaminants  
                ).reads

                // Combine paired and single-ends into a single channel
                ch_shortreads_preprocessed = ch_shortreads_paired_preprocessed.mix(ch_shortreads_single_preprocessed)
                
                // Collect versions from BBDUK
                ch_versions = ch_versions.mix( BBMAP_BBDUK_PAIRED.out.versions )
                ch_versions = ch_versions.mix( BBMAP_BBDUK_SINGLE.out.versions )

            } else {
                // No BBDUK preprocessing, just pass through input reads to FASTQC
                ch_shortreads_preprocessed = ch_verified_samplesheet.paired.mix(ch_verified_samplesheet.single)
            }

            // If the qc_tool is 'fastp', run fastp for QC
            if (tool == 'fastp') {
                
            }

            //
            // MODULE: Run BBmerge
            //
            // Initialise channel
            final_input_reads = Channel.empty()

            if (params.bbmerge_pairs ) {

                BBMAP_BBMERGE_PRE ( BBMAP_BBDUK_PAIRED.out.reads, [] )
                .merged
                    .map {
                        sample, reads ->
                            return [ sample + [single:true ], reads.flatten() ] // Flatten the reads list for each sample
                    }
                    .set { ch_merged_reads_pe }
                ch_versions = ch_versions.mix(BBMAP_BBMERGE_PRE.out.versions.first())


                //ch_merged_reads = ch_merged_reads_pe.mix(ch_shortreads_single_preprocessed) // unneeded as else statement takes care if not paired-end or merging not selected

                // Same error for code below
                // ch_merged_reads = BBMAP_BBMERGE_PRE (ch_shortreads_preprocessed, []).merged   

                
                // If not merging runs, use bbmerged reads as final input for metaphlan
                if ( !params.perform_runmerging ) {
                    final_input_reads = ch_merged_reads_pe 
                }
            } else {
                
            // Ensure `final_input_reads` is assigned even if `bbmerge_pairs` is false
            final_input_reads = ch_shortreads_preprocessed
            }

            /*
                Run host removal, run_merge, etc... to be completed
            */

            //
            // MODULE: Run FastQC post-processing
            //
            if ( tool == 'fastqc' ) {    
                FASTQC_POST ( 
                    ch_shortreads_preprocessed 
                )
                // Mix any output files for MultiQC into a multiqc channel
                ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.collect{it[1]})
                ch_versions = ch_versions.mix(FASTQC_POST.out.versions.first())
            }
        }
    }
    
    emit:
    final_input_reads
    ch_multiqc_files
    versions = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
