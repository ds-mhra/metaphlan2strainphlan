<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-metaphlanstrainphlan_logo_dark.png">
    <img alt="nf-core/metaphlan2strainphlan" src="docs/images/nf-core-metaphlanstrainphlan_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/metaphlanstrainphlan/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/metaphlanstrainphlan/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/metaphlanstrainphlan/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/metaphlanstrainphlan/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/metaphlan2strainphlan)


## Introduction

**nf-core/metaphlan2strainphlan** is a bioinformatics pipeline that ...


This pipeline runs both MetaPhlAn 4.1 and StrainPhlAn 4.0 from the biobakery collection to analyse metagenomic shotgun sequencing data. The former allows users to profile the composition of microbial communities while the latter characterises sample sets to a strain-level resolution of the species of interest.

Please refer to the tools' Github repository:
- https://github.com/biobakery/biobakery/wiki/MetaPhlAn-4.1
- https://github.com/biobakery/biobakery/wiki/strainphlan4


<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->

## Default steps
1. Pre QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Optional pre-processing of FASTQ using ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | [`fastp`](https://github.com/OpenGene/fastp) | [`BBDuk`](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)) which includes:
    - trimming and quality filtering 
    - adapter or contamination removal
    - merging any paired reads 
3. Post QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Download database and perform taxonomic classification and profiling with ([`MetaPhlAn4.1`](https://github.com/biobakery/biobakery/wiki/MetaPhlAn-4.1))
5. Merge MetaPhlAn results and standardise into an output table (`METAPHLAN_MERGEMETAPHLANTABLES`)
6. Generate abundance tables of species
7. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

Additional steps:
8. Characterise strains using ([`StrainPhlAn4`](https://github.com/biobakery/biobakery/wiki/strainphlan4))
 

## Using Nextflow

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Before running the pipeline, please refer to the [usage documentation](docs/usage.md) for more details and further functionality.

Now, run using:

```bash
nextflow run nf-core-metaphlanstrainphlan \
	-profile docker,gcb \
  --samplesheet samplesheet.csv \
	--installdb --metaphlan_index mpa_vJun23_CHOCOPhlAnSGB_202403 \
	--run_metaphlan --run_strainphlan \
	--bbmerge_pairs \
	--qc_tool 'bbduk,fastqc' \
  --outdir './results'
```


This pipeline can be launched with different profiles such as `docker` and `gcb` configurations. 

> [!WARNING]
> **NOTE: the `gcb` profile in the nextflow.config will need to be amended to match user's own google credentials. E.g. `bucket_name` and `project_id` in the file should be filled in.**

See below for more information about profiles.


> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).
<br>

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/metaphlan2strainphlan/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/metaphlan2strainphlan/output).

## Credits

nf-core/metaphlanstrainphlan was originally written by Dammy Shittu.


## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/metaphlan2strainphlan for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the nf-core publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
