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

include { PREPROCESSING                      } from './workflows/preprocessing'
include { PROFILING; METAPHLAN_TO_STRAINPHLAN; STRAIN_CHARACTERISATION } from './workflows/metaphlanstrainphlan'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_metaphlanstrainphlan_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_metaphlanstrainphlan_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_metaphlanstrainphlan_pipeline'


//
// Print help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run nf-core-metaphlanstrainphlan <ARGUMENTS>

Arguments:

  Input:
  --samplesheet             csv file containing .fastq or .fastq.gz files with header sample,fastq_1,fastq_2
  --input_folder........            Directory containing file pairs: may need to override the following parameters: file_spacer,file_suffix

  Preprocess sample files:
  --perform_shortread_qc    Performs shortread quality control, not necessary if `--qc_tool` is already used.

  --qc_tool                 A comma-separated list for quality control review, must be used with `--perform_shortread_qc`.
                            options: bbduk|fastqc|fastp       example: --qc_tool 'bbduk,fastqc'
  
  --shortread_qc_contaminantslist       fa|fasta|csv|tsv|txt
  
  --bbmerge_pairs           Merge paired-end fastq files using BBMerge.

  MetaPhlAn:
  --installdb       
  
  --run_metaphlan           Run MetaPhlAn.
  
  --metaphlan_index         Name of available databases ie `mpa_vJun23_CHOCOPhlAnSGB_202403`.
  
  --metaphlan_db            Path to reference database, must not be used if installing.
                            (must include the prefix shared across reference files)

  StrainPhlAn:
  --run_strainphlan         Run StrainPhlAn

  --strainphlan_clades 
  
  --reference_genomes       csv file containing a list of all reference genomes/ strains to be used with StrainPhlAn
                            The headers must be: ref_genome_strainID,ref_genome_path,species_name,taxonID,associated_clade, where the first two mandatory
                            example: GCA_008121495.1,path/to//GCA_008121495.1_ASM812149v1_genomic.fna,Ruminococcus_gnavus,33038,t__SGB4584

  --metadata                txt file

  --strainphlan_db          Optional, path. Running MetaPhlAn and excluding this parameter allows StrainPhlAn to use the same downloaded database.

  --clade_markers           Optional, path. Running MetaPhlAn and excluding this parameter allows StrainPhlAn to use the output of MetaPhlAn.

  Output Location:
  --outdir                  Output folder 

  Skip processes:
  --skip_preprocessing_qc


""".stripIndent()
}


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
    PREPROCESSING (
        ch_samplesheet,
    )

    // Run profiling with MetaPhlAn
    PROFILING (
        PREPROCESSING.out.final_input_reads,
    )

    // Extract clades from metaphlan results; download genomes w NCBI datasets
    METAPHLAN_TO_STRAINPHLAN (
        PROFILING.out.ch_profiles
    )

    // Run characterisation of strains via StrainPhlAn
    STRAIN_CHARACTERISATION (
        PROFILING.out.ch_sam_files, 
        PROFILING.out.ch_final_dbs,
        METAPHLAN_TO_STRAINPHLAN.out.clades_list, 
        []
        // all_references
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
