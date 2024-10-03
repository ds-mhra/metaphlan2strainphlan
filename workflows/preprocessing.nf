/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UNTAR                  } from '../modules/nf-core/untar/main'
//include { INPUT_CHECK            } from '../subworkflows/input_check'
include { FASTP as FASTP_PAIRED  } from '../modules/nf-core/fastp/main'
include { FASTP as FASTP_SINGLE  } from '../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_PRE   } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_POST  } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_metaphlanstrainphlan_pipeline'
include { BBMAP_ALIGN            } from '../modules/nf-core/bbmap/align/main' 
include { BBMAP_BBDUK as BBMAP_BBDUK_PAIRED; BBMAP_BBDUK as BBMAP_BBDUK_SINGLE } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBMERGE          } from '../modules/nf-core/bbmap/bbmerge/main' 


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
    WORKFLOW: PREPROCESSING
*/

// Initialise empty channels


workflow PREPROCESSING {

    take:
    ch_samplesheet     // channel: samplesheet read in from --input 
    
    main:
    ch_versions                   = Channel.empty()
    ch_multiqc_files              = Channel.empty()
    ch_shortreads_pe_preprocessed = Channel.empty()
    ch_shortreads_se_preprocessed = Channel.empty()

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
                // Map sample data with file checks by creating a tuple of sample metadata and fastq_files
                def fastq_files = [
                    row.fastq_1 ? file(row.fastq_1, checkIfExists: true) : error("Missing fastq_1 file for sample ${row.sample}"),    // Check for, and always include, fastq_1 (this includes SE and PE samples) 
                ]

                // Ensure row has an 'id' field
                def id = row.sample 

                // Check if fastq_2 is available; if not, handle it accordingly
                if (row.fastq_2) {
                    fastq_files << file(row.fastq_2, checkIfExists: true)           // Add fastq_2 if present
                    return [id: row.sample, fastq_files: fastq_files, type: "paired"]                      // Return with a flag for paired reads
                } else {
                    return [id: row.sample, fastq_files: fastq_files, type: "single"]                      // Return with a flag for single reads
                }

                // Return tuple ????????????????
                //tuple(row.sample, fastq_files)
            }

        // Separate the channels based on the read type using a conditional structure
        ch_samplesheet_single = ch_verified_samplesheet
            .filter { it[2] == "single" }           // Filter and emit only single-end reads
            .map { [it[0], it[1][0]] }              // Keep only sample and fastq_1

        ch_samplesheet_paired = ch_verified_samplesheet
            .filter { it[2] == "paired" }           // Filter and emit paired-end reads
            .map { [it[0], it[1][0], it[1][1]] }    // Keep sample, fastq_1, and fastq_2

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
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PRE.out.zip.collect{it[1]})             

    // Perform QC if either perform_shortread_qc or a qc_tool is selected
    if (params.perform_shortread_qc || params.qc_tool) {

        // Error message
        // Define the allowed tools list and esure it is treated as a list
        def allowed_qc_tools = ['bbduk', 'fastp', 'fastqc', 'bbmerge']
        def selected_qc_tools = params.qc_tool instanceof List ? params.qc_tool : [params.qc_tool]

        // After treating each tool like a list, validate that all selected tools are supported
        selected_qc_tools.each { tool ->
            if (!allowed_qc_tools.contains(tool)) {
                error "[metaphlan2strainphlanPipeline] ERROR: Unsupported QC tool has been selected: ${tool}"
            }
        }

        // If the qc_tool is 'bbduk', run BBDUK for QC trimming and filtering
        if (params.qc_tool.contains('bbduk')) {

            contaminants = params.shortread_qc_contaminantslist ? file(params.shortread_qc_contaminantslist) : []

            if (params.shortread_qc_contaminantslist && !contaminants.extension.matches(".*(csv|tsv|txt)")) {
                error "[metaphlan2strainphlanPipeline] ERROR: Contaminants or adapter list requires a different format and/or extension. Check input: --shortread_qc_contaminantslist ${params.shortread_qc_contaminantslist}"
            }

            //
            // MODULE: Run BBMAP_BBDUK
            // Trim and filter paired-end reads with BBDUK
            ch_shortreads_pe_preprocessed = BBMAP_BBDUK_PAIRED( 
                ch_samplesheet_paired, contaminants
            ).reads
            
            // Trim and filter single-end reads with BBDUK
            ch_shortreads_se_preprocessed = BBMAP_BBDUK_SINGLE (
                ch_samplesheet_single, contaminants  
            ).reads

            // Combine paired and single-ends into a single channel
            ch_shortreads_preprocessed = ch_shortreads_pe_preprocessed.mix(ch_shortreads_se_preprocessed)
            
            // Collect versions from BBDUK
            ch_versions = ch_versions.mix( BBMAP_BBDUK_PAIRED.out.versions )
            ch_versions = ch_versions.mix( BBMAP_BBDUK_SINGLE.out.versions )

        } else {
            // No BBDUK preprocessing, just pass through input reads to FASTQC
            ch_shortreads_preprocessed = ch_samplesheet_paired.mix(ch_samplesheet_single)
        }

        //
        // MODULE: Run BBmerge for paired-end reads
        //
        // Initialise channel
        final_input_reads = Channel.empty()

        if (params.bbmerge_pairs) {

            BBMAP_BBMERGE ( ch_shortreads_pe_preprocessed, [] ).merged
                // .map {
                //     sample, reads ->
                //         return [ sample + [single:true ], reads.flatten() ] // Flatten the reads list for each sample
                // }
                .set { ch_merged_reads_pe }            // ch_merged_reads = BBMAP_BBMERGE (ch_shortreads_preprocessed, []).merged   

            //ch_merged_reads = ch_merged_reads_pe.mix(ch_shortreads_single_preprocessed) // unneeded as else statement takes care if not paired-end or merging not selected

            // Mix any output files 
            ch_versions = ch_versions.mix(BBMAP_BBMERGE.out.versions.first())
            
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


        // If the qc_tool is 'fastp', run fastp for QC trimming, filtering and merging
        if (params.qc_tool.contains('fastp')) {
            ch_fastp_pe = FASTP_PAIRED ( 
                ch_samplesheet_paired, contaminants 
            )
            ch_fastp_se = FASTP_SINGLE ( 
                ch_samplesheet_single, contaminants 
            )
        // Mix any output files
        ch_multiqc_files = ch_multiqc_files
            .mix(FASTP_PAIRED.out.zip.collect{it[1]})
            .mix(FASTP_SINGLE.out.zip.collect{it[1]})
        
        ch_versions = ch_versions
            .mix(FASTP_PAIRED.out.versions.first())
            .mix(FASTP_SINGLE.out.versions.first())
        
        final_input_reads = ch_fastp_pe.reads.mix(ch_fastp_se.reads)
        }

        //
        // MODULE: Run FastQC post-processing
        //
        if (params.qc_tool.contains('fastqc')) { 
            FASTQC_POST ( 
                final_input_reads 
            )
            // Mix any output files for MultiQC into a multiqc channel
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(FASTQC_POST.out.versions.first())
        }
    } else {
        final_input_reads = ch_verified_samplesheet.paired.mix(ch_verified_samplesheet.single)
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
