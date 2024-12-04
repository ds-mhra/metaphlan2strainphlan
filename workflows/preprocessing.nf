//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include { INPUT_CHECK            } from '../subworkflows/input_check'
include { UNTAR                  } from '../modules/nf-core/untar/main'
include { FASTP as FASTP_PAIRED; FASTP as FASTP_SINGLE  } from '../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_PRE; FASTQC as FASTQC_POST   } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_metaphlanstrainphlan_pipeline'
include { BBMAP_ALIGN            } from '../modules/nf-core/bbmap/align/main' 
include { BBMAP_BBDUK as BBMAP_BBDUK_PAIRED } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBMAP_BBDUK_SINGLE } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBMERGE          } from '../modules/nf-core/bbmap/bbmerge/main' 


/*
    WORKFLOW: PREPROCESSING
*/

// Initialise empty channels


workflow PREPROCESSING {

    take:
    ch_samplesheet     // ch_samplesheet or ${params.samplesheet}; read in from --input or --samplesheet
    
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
                params.samplesheet,
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
                // Map sample data with file checks by creating a tuple of sample metadata (id) and fastq_files
                def id = row.sample 

                def fastq_files = [
                    row.fastq_1 ? row.fastq_1 : error("Missing fastq file for sample ${id}"),    // Check for, and always include, fastq_1 (this includes SE and PE samples) 
                ]

                // Check if fastq_2 is available; if not, handle it accordingly and return list with 3 elements, 'sample', 'fastq, 'type'
                if (row.fastq_2) {
                    fastq_files << row.fastq_2                                  // Append fastq_2 if present
                }
                // Debugging: Print the structure before returning
                println("Sample ID: ${id}, Fastq Files: ${fastq_files}")
                return [id: id, fastq_files: fastq_files]       
            }
        
        // Modified channel splitting
        ch_samplesheet_paired = ch_verified_samplesheet
            .filter { it.fastq_files?.size() == 2 }                         // ensures that fastq_2 exists 
            .map { [[ id: it.id, single_end: false ], it.fastq_files] }     // since modules use {meta.id} to access the ID, the ID must be wrapped (nested) in the meta map while adding a single flag
            .ifEmpty { log.warn "No paired-end samples found." }

        ch_samplesheet_single = ch_verified_samplesheet
            .filter { it.fastq_files?.size() == 1 }                         // ensures only fastq_1 exists
            .map { [[ id: it.id, single_end: true ], it.fastq_files] }      
            .ifEmpty { log.warn "No single-end samples found." }
<<<<<<< HEAD
        
/*
    - indices don't work??try??   .filter { it[1].size() }         OR      .map { [[it[0], single_end: true], it[1]] }      AND      .map { [[it[0], single_end: false], it[1]] }
*/


/*  ~ ORIGINAL, not working, wont properly call indices ~
        // Separate the channels based on the read type using a conditional structure
        ch_samplesheet_single = ch_verified_samplesheet
            .filter { it[2] == "single" }           // Filters by second (0,1,2) element in the tuple, ie type, to emit single-end reads
            .map { [it[0], it[1][0]] }              // Keep only sample and fastq_1

        ch_samplesheet_paired = ch_verified_samplesheet
            .filter { it[2] == "paired" }           // Filter by second (0,1,2) element array in tuple, type, to emit paired-end reads
            .map { [it[0], it[1][0], it[1][1]] }    // Keeps sample, fastq_1, and fastq_2 
*/        
=======
>>>>>>> preprocess_metaphlan
    
            // Debugging channel structure
            // ch_verified_samplesheet.view { "Original: $it" }
            // ch_samplesheet_paired.view { "Debug BBDUK PE input: meta=${it[0]}, reads=${it[1]}" }
            // ch_samplesheet_single.view { "Debug BBDUK SE input: meta=${it[0]}, reads=${it[1]}" }


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
    // Filter out null or empty tuples in the paired-end and single-end channels, combine and use as input to fastqc
    ch_samplesheet_paired = ch_samplesheet_paired.filter { it != null }
    ch_samplesheet_single = ch_samplesheet_single.filter { it != null }

    FASTQC_PRE ( 
        ch_samplesheet_paired.mix(ch_samplesheet_single)
    )
    ch_versions = ch_versions.mix(FASTQC_PRE.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PRE.out.zip.collect{it[1]})             


    // Perform QC if either perform_shortread_qc or a qc_tool is selected
    if (params.perform_shortread_qc || params.qc_tool) {

        // Error message:
        // Define the allowed tools list and esure it is treated as a list
        def allowed_qc_tools = ['bbduk', 'fastp', 'fastqc', 'bbmerge']
    
        // Convert params.qc_tool to a list, and handles different input formats on CLI
        def selected_qc_tools = []
        if (params.qc_tool instanceof String) {
            // Handle comma-separated string and splits it up
            selected_qc_tools = params.qc_tool.split(',').collect { it.trim() }
        } else if (params.qc_tool instanceof List) {
            // Handle list as input
            selected_qc_tools = params.qc_tool
        } else if (params.qc_tool != null) {
            selected_qc_tools = [params.qc_tool.toString()]
        }

        // After treating each tool like a list, validate that all selected tools are supported
        selected_qc_tools.each { tool ->
            if (!allowed_qc_tools.contains(tool)) {
                error "[metaphlan2strainphlanPipeline] ERROR: Unsupported QC tool has been selected: ${tool}"
            }
        }

        // If the qc_tool is 'bbduk', run BBDUK for QC trimming and filtering
        if (params.qc_tool.contains('bbduk')) {
            
            println "Processing with BBDUK..."
            println "Debug: QC tool selection: ${selected_qc_tools}"    // shows list
<<<<<<< HEAD
                // println "Debug: QC tool selection: ${params.qc_tool}"       // shows string
                // println "Debug: Contains bbduk?  ${selected_qc_tools.contains('bbduk')}"
                // println "Debug: Contains fastp?  ${selected_qc_tools.contains('fastp')}"
                // println "Debug: Contains fastqc? ${selected_qc_tools.contains('fastqc')}"
            
/*          Processing with BBDUK...
            Sample ID: Day_1A, Fastq Files: [gs://bucket_name/y/Day1_1A_S38_R1_001.fastq.gz, gs://bucket_name/y/Day1_1A_S38_R2_001.fastq.gz]
            Sample ID: Day_1B, Fastq Files: [gs://bucket_name/y/Day1_1B_S3_R1_001.fastq.gz]
            Debug: QC tool selection: bbduk,fastqc
            Debug: QC tool selection: [bbduk, fastqc]
            Debug: Contains bbduk? true
            Debug: Contains fastp? false
            Debug: Contains fastqc? true
*/
=======
   
>>>>>>> preprocess_metaphlan

            // If contaminants file is absent, use an empty list
            contaminants = params.shortread_qc_contaminantslist ? file(params.shortread_qc_contaminantslist) : []

            if (params.shortread_qc_contaminantslist && !contaminants.extension.matches(".*(fa|fasta|csv|tsv|txt)")) {
                error "[metaphlan2strainphlanPipeline] ERROR: Contaminants or adapter list requires a different format and/or extension. Check input: --shortread_qc_contaminantslist ${params.shortread_qc_contaminantslist}"
            }

            // Since Nextflow keeps interpreting and removing the gs://bucket_name in the file paths in favour of local directories
            // Store absolute paths before BBDUK by converting to string and defining bucket_paths
            ch_samplesheet_paired_tracked = ch_samplesheet_paired.map { id, fastq_files ->
<<<<<<< HEAD
                id.original_paths = fastq_files.collect { it.toString() }            // modify id object by adding original_paths,      OR id.bucket_paths OR def bucket_paths 
                [id, fastq_files]
            }
                // DEBUG works fine - may need to change id.bucket_paths to id and/or fastq_files to bucket_paths
                // ch_samplesheet_paired_tracked.view { id, fastq_files ->
                //     println "Debug 1 PE - PAIRED Check sample: ${id.id}, and its bucket paths: ${id.original_paths}"
                //     println "Debug 1 PE - PAIRED Current files: ${fastq_files}"
                // }
                // ch_samplesheet_paired_tracked.view { 
                //     println "Debug 1.1  PE - PAIRED; ALL ${it}"
                // }
=======
                id.original_paths = fastq_files.collect { it.toString() }            // modify id object by adding original_paths 
                [id, fastq_files]
            }
                
>>>>>>> preprocess_metaphlan

            //single
            ch_samplesheet_single_tracked = ch_samplesheet_single.map { id, fastq_file ->
                id.original_path = fastq_file.collect { it.toString() }
                [id, fastq_file]
            }
<<<<<<< HEAD
                // ch_samplesheet_single_tracked.view { id, fastq_file ->
                //     println "Debug 1 SE - Check SINGLE sample: ${id.id}, and its bucket path: ${id.original_path}"
                //     println "Debug 1 SE - Current file: ${fastq_file}"
                // }
                // ch_samplesheet_single_tracked.view { 
                //     println "Debug 1.1 SE - SINGLE; ALL ${it}"
                // }
=======
             
>>>>>>> preprocess_metaphlan

            // MODULE: Run BBMAP_BBDUK
            // Trim and filter paired- and single-end reads with BBDUK
            BBMAP_BBDUK_PAIRED( 
                ch_samplesheet_paired_tracked, contaminants
            )
            ch_shortreads_pe_preprocessed = BBMAP_BBDUK_PAIRED.out.reads
<<<<<<< HEAD
            // ch_shortreads_pe_preprocessed.view { "BBDUK PE output: $it" }
=======
            ch_shortreads_pe_preprocessed.view { "BBDUK PE output: $it" }
>>>>>>> preprocess_metaphlan

            BBMAP_BBDUK_SINGLE( 
                ch_samplesheet_single_tracked, contaminants
            )
            ch_shortreads_se_preprocessed = BBMAP_BBDUK_SINGLE.out.reads
<<<<<<< HEAD
            // ch_shortreads_se_preprocessed.view { "BBDUK SE output: ${it}" } 
                // ch_shortreads_se_preprocessed.view { id, bucket_path ->
                //     println "Debug 2 - SINGLE Check sample: ${id.id}, and its bucket path: ${bucket_path}"
                //     println "Debug 2 - SINGLE Current file: ${bucket_path}"
                // }
=======
            ch_shortreads_se_preprocessed.view { "BBDUK SE output: ${it}" } 
          
>>>>>>> preprocess_metaphlan
            // Combine paired and single-ends into a single channel
            ch_shortreads_preprocessed = ch_shortreads_pe_preprocessed.mix(ch_shortreads_se_preprocessed)
            
            // Collect versions from BBDUK
            ch_versions = ch_versions.mix( BBMAP_BBDUK_PAIRED.out.versions )
            ch_versions = ch_versions.mix( BBMAP_BBDUK_SINGLE.out.versions )
<<<<<<< HEAD
/*
            //paired
            Debug 1 PE - PAIRED Check sample: Day_1A, and bucket paths: [gs://bucket_name/y/Day1_1A_S38_R1_001.fastq.gz, gs://bucket_name/y/Day1_1A_S38_R2_001.fastq.gz]
            Debug 1 PE - PAIRED Current files: [gs://bucket_name/y/Day1_1A_S38_R1_001.fastq.gz, gs://bucket_name/y/Day1_1A_S38_R2_001.fastq.gz]
            Debug 1.1 - PAIRED; ALL [[id:Day_1A, single_end:false, bucket_paths:[gs://bucket_name/y/Day1_1A_S38_R1_001.fastq.gz, gs://bucket_name/y/Day1_1A_S38_R2_001.fastq.gz]], [gs://bucket_name/y/Day1_1A_S38_R1_001.fastq.gz, gs://bucket_name/y/Day1_1A_S38_R2_001.fastq.gz]]
                // using:  it[1]    [gs://bucket_name/y/Day1_1A_S38_R1_001.fastq.gz, gs://bucket_name/y/Day1_1A_S38_R2_001.fastq.gz]
                // using:  it[0]    [id:Day_1A, single_end:false, bucket_paths:[gs://bucket_name/y/Day1_1A_S38_R1_001.fastq.gz, gs://bucket_name/y/Day1_1A_S38_R2_001.fastq.gz]]
            BBDUK PE output: [[id:Day_1A, single_end:false, bucket_paths:[gs://bucket_name/y/Day1_1A_S38_R1_001.fastq.gz, gs://mhra-ngs-dev-ut8t-training/lekshmi_bucket/152_Day1_1A_S38_R2_001.fastq.gz]], [/workdir/3f/afd9bdae6a5692c2784f123297f1f1/Day_1A_1.fastq.gz, /workdir/3f/afd9bdae6a5692c2784f123297f1f1/Day_1A_2.fastq.gz]]

            //single
            Debug 1 SE - Check SINGLE sample: Day_1B, and its bucket path: [gs://bucket_name/y/Day1_1B_S3_R1_001.fastq.gz]
            Debug 1 SE - Current file: [gs://bucket_name/y/Day1_1B_S3_R1_001.fastq.gz]
            Debug 1.1 SE - SINGLE; ALL [[id:Day_1B, single_end:true, bucket_path:[gs://bucket_name/y/Day1_1B_S3_R1_001.fastq.gz]], [gs://bucket_name/y/Day1_1B_S3_R1_001.fastq.gz]]
                // using:  it[1]     [gs://bucket_name/y/Day1_1B_S3_R1_001.fastq.gz]
                // using:  it[0]      [id:Day_1B, single_end:true, bucket_path:[gs://bucket_name/y/Day1_1B_S3_R1_001.fastq.gz]]
            BBDUK SE output: [[id:Day_1B, single_end:true, bucket_path:[gs://bucket_name/y/Day1_1B_S3_R1_001.fastq.gz]], /workdir/f4/63d84f7ec04c5b84305cd770a41db6/Day_1B.fastq.gz]
                // using:  it[0]      [id:Day_1B, single_end:true, bucket_path:[gs://bucket_name/y/Day1_1B_S3_R1_001.fastq.gz]]

            //all
            Debug: FINAL final_input_reads set: [[id:Day_1B, single_end:true, bucket_path:[gs://bucket_name/y/Day1_1B_S3_R1_001.fastq.gz]], /workdir/f4/63d84f7ec04c5b84305cd770a41db6/Day_1B.fastq.gz]
            ch_metaphlan_input hereeeeee [[id:Day_1B, single_end:true, bucket_path:[gs://bucket_name/y/Day1_1B_S3_R1_001.fastq.gz]], /workdir/f4/63d84f7ec04c5b84305cd770a41db6/Day_1B.fastq.gz]
*/
=======

>>>>>>> preprocess_metaphlan

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
                .map { id, merged_fastq ->                              // change PE label on merged sample using `.map { }` so that it is marked as paired_end and merged
<<<<<<< HEAD
                    [id + [merged: true], merged_fastq]                 // OR merged_fastq.flatten()    // Flatten reads list for each sample, should remove null sample after merging
                }
                .set { ch_merged_reads_pe }
            ch_merged_reads_pe.view { "Debug: ch_merged_reads_pe from BBMERGE set: $it" }       // [[id:Day_1A, single_end:false,  merged:true], /workdir/4e/d38889cca9c2dc16193a4d0e0a5b84/Day_1A_merged.fastq.gz]
=======
                    [id + [merged: true], merged_fastq]                 
                }
                .set { ch_merged_reads_pe }
            ch_merged_reads_pe.view { "Debug: ch_merged_reads_pe from BBMERGE set: $it" }   
>>>>>>> preprocess_metaphlan

            final_input_reads = ch_merged_reads_pe.mix(ch_shortreads_se_preprocessed)
                
            // Mix any output files 
            ch_versions = ch_versions.mix(BBMAP_BBMERGE.out.versions.first())
            
            // If not merging runs, use bbmerged reads as final input for metaphlan
            if ( !params.perform_runmerging ) {
                final_input_reads = final_input_reads 
            } else {
                // ... perform run merging...
            }
        } else {
            
        // Ensure `final_input_reads` is assigned even if `bbmerge_pairs` is false
        final_input_reads = ch_shortreads_preprocessed
        }
<<<<<<< HEAD
        // final_input_reads.view { "Debug: FINAL final_input_reads set: $it" }      // [[id:Day_1A, single_end:false, merged:true], /workdir/4e/d38889cca9c2dc16193a4d0e0a5b84/Day_1A_merged.fastq.gz]
        //                                                                             [[id:Day_1B, single_end:true], /workdir/84/1cd89ef46c007df4b41c8fe12a5711/Day_1B.fastq.gz]

=======
        final_input_reads.view { "Debug: FINAL final_input_reads set: $it" } 
>>>>>>> preprocess_metaphlan

        /*
            Run host removal, run_merge, etc... to be completed
        */


        // If the qc_tool is 'fastp', run fastp for QC trimming, filtering and merging
        // NB. fastp shouldn't be run w bbmerge
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
        } else {
            // No fastp preprocessing, just pass through input reads
            ch_shortreads_preprocessed = ch_samplesheet_paired.mix(ch_samplesheet_single)
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
<<<<<<< HEAD
        final_input_reads = ch_samplesheet_paired.mix(ch_samplesheet_single)        // changed from ch_verified_samplesheet.paired.mix(ch_verified_samplesheet.single)
=======
        final_input_reads = ch_samplesheet_paired.mix(ch_samplesheet_single)       
>>>>>>> preprocess_metaphlan
    }
    
    emit:
    final_input_reads
    ch_multiqc_files
    versions = ch_versions 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
