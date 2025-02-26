# nf-core/metaphlanstrainphlan: Usage

## :warning: Please read this documentation before running the pipeline.

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This pipeline runs both MetaPhlAn 4.1 and StrainPhlAn 4.0 from the biobakery collection to analyse metagenomic shotgun sequencing data. The former allows users to profile the composition of microbial communities while the latter characterises sample sets to a strain-level resolution of the species of interest.

Please refer to the tools' Github repository:
- https://github.com/biobakery/biobakery/wiki/MetaPhlAn-4.1
- https://github.com/biobakery/biobakery/wiki/strainphlan4


## Samplesheet input

Create a samplesheet with information about the samples to be analysed prior to running the pipeline. This can be done manually or by running the *metaphlan_samplesheet.sh* script. 

```bash
bash nf-core-metaphlanstrainphlan/bin/metaphlan_samplesheet.sh "[fastq directory]" "[samplesheet name]"
```
<br>

The final samplesheet file must be a comma-separated (.csv) file with at least 3 columns, and a header row as shown in the example below. 

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet, simply leave `fastq_2` empty.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
SAMPLE1,SAMPLE1_S1_L002_R1_001.fastq.gz,SAMPLE1_S1_L002_R2_001.fastq.gz
SAMPLE2,SAMPLE2_S1_L003_R1_001.fastq.gz,SAMPLE2_S1_L003_R2_001.fastq.gz
SAMPLE3,SAMPLE3_S1_L004_R1_001.fastq.gz,
```
<br>

Use the `--samplesheet` parameter to specify its location. 

```bash
--samplesheet '[path to samplesheet file]'
```
<br>

### Multiple runs of the same sample

The `sample` identifiers have to be the same when users have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```


| Column    |Description             |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Ensure any spaces in sample names are converted to underscores (_). |
| `fastq_1` | Full path to FastQ file for reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
<br>

## Genomes list input
If `--run_strainphlan` is used, genome reference files can be included in the StrainPhlAn analysis.
The optional genome_list file must be a comma-separated (.csv) file with at least 5 columns, and a header row as shown in the example below. 

The file contains relevant metadata in regards to the reference species/strains with and can be done manually or by running the *list_genomes.sh* script. The final two fields must be filled manually.

NB. 
- All fields are compulsory, except the `taxonID`. 
- This pipeline specifically looks for *.fna* files in a directory structure similar to that of NCBI downloads (please see below).

```bash
bash nf-core-metaphlanstrainphlan/bin/list_genomes.sh '[path to reference_1]' '[path to reference_2]' '[path to genome_list file].csv'
```
<br>

```csv title="genome_list.csv"
ref_genome_strainID,ref_genome_path,species_name,taxonID,associated_clade
GCA_000123456.1,Escherichia coli_ncbi_dataset/data/GCA_000123456.1/GCA_000123456.1_ASM12345v1_genomic.fna,Escherichia coli,,t__SGB10068
```


| Column    |Description             |
| --------- | ----------------------------------------------------- |
| `ref_genome_strainID`  | Strain ID name. Can be found in the file name. This entry will be unique for each strain listed. |
| `ref_genome_path` | Full path to the strain's reference genome ".fna". This can be local or in cloud storage. |
| `species_name` | Provide full scientific name.   |
| `taxonID` | NCBI Taxonomy ID number.  |
| `associated_clade` | Associated clade name, can be found in MetaPhlAn results or a taxonomy database.  |

### Reference genome directory structure
Based on NCBI downloads.
```
genome_directory/
  └── species_name_ncbi_dataset/
      └── data/
          └── strain_name/
              └── strain_name.fna
```
<br>

## Running the pipeline

The typical default command for running the pipeline is as follows:

```bash
nextflow run nf-core-metaphlanstrainphlan \
	-profile <docker/singularity/gcb/...> \
  --samplesheet samplesheet.csv \
	--installdb --metaphlan_index mpa_vJun23_CHOCOPhlAnSGB_202403 \
	--run_metaphlan --run_strainphlan \
	--bbmerge_pairs \
	--qc_tool 'bbduk,fastqc' \
  --outdir './results'
```

This will launch the pipeline with the `docker` and `gcb` configuration profile. 

**NOTE: the `gcb` profile in the nextflow.config will need to be amended to match user's own google credentials. E.g. `bucket_name` and `project_id` in the file should be filled in.**

See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Preprocessing:
QC_TOOLS available: `fastqc`, `fastp`, `bbduk`

The pipeline is able to intake single- and paired-end reads, pre-process the files according to the selected quality control (QC) tool, and subsequently merge any pairs with BBMap, via `--bbmerge_pairs`.

### Running MetaPhlAn and StrainPhlAn

If `--installdb` is selected, users will need to specify the index database name, ie `--metaphlan_index mpa_vJun23_CHOCOPhlAnSGB_202403`, if this is not specified, the latest version is downloaded automatically.

If an index database is already present, simply use `--metaphlan_db` to direct to the top-level directory containing the index. The pipeline will search for the *.PKL* file within this directory, however, if `--strainphlan_db` is specified, it must point to the *.PKL* file within this directory.

Note that MetaPhlAn can be run without StrainPhlAn. The preparation for StrainPhlAn (e.g. creating consensus markers and extracting clade markers) can be skipped with `--skip_strainphlan_prep` however, the directory containing these files must be specified in addition to the clade of interest via `--strainphlan_clades` and the reference genome list.


## Parameters:
```
usage: nextflow run nf-core-metaphlanstrainphlan
	[--bbmerge_pairs]
	[--qc_tool QC_TOOLS] [--shortread_qc_contaminantslist ADAPTERS_FILE] [--perform_shortread_qc] [--skip_preprocessing_qc]
  [--run_metaphlan] [--run_strainphlan] [--skip_strainphlan_prep]
	[--installdb] [--metaphlan_index MPA_INDEX_NAME]
  [--metaphlan_db DATABASE_DIR] [--strainphlan_db DATABASE_PKL]
  [--reference_genomes GENOME_LIST_FILE]
  [--strainphlan_clades CLADE_NAME]
  [--consensus_markers MARKER_DIR] [--clade_markers CLADE_DIR] [--metadata STRAINPHLAN_METADATA_FILE]
	[--graphlan_plots]
  [-profile <docker,googlebatch>]
  [-h]
  [--samplesheet INPUT_FILE] [OUTPUT_DIR]

```
<br>

*Note [graphlan_annotate.py](../bin/graphlan_annotate.py) and [graphlan.py](../bin/graphlan.py) scripts can be used.*

<br>

### -params-file
To repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, users can supply a params file in a `yaml` or `json` format via `-params-file <file>`.

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/metaphlanstrainphlan -profile docker -params-file params.yaml
```

with `params.yaml` containing something like:

```yaml
samplesheet: './samplesheet.csv'
outdir: './results/'
bbmerge_pairs: true
<...>
```

:warning:
Do not use `-c <file>` to specify parameters as this may result in errors. Custom config files specified with `-c` should only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
<br>



## Extra information
### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/metaphlanstrainphlan
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data to ensure a specific version of the pipeline code and software are used.

First, go to the Github page releases directory and find the latest pipeline version - numeric only (eg. `1.3.1`) if available. Specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1` and it will be logged in the output reports.


Tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

### Core Nextflow arguments

#### `-profile`

Use this parameter to choose a configuration profile.

:::info
Using Docker or Singularity containers is highly recommend for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

#### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.


#### Running in the background

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.


#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
