# nf-ncov-voc

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)


## Introduction

**nf-ncov-voc** is a bioinformatics analysis workflow used to perform variant calling for SARS-CoV-2 genomes to identify and profile mutations in Variants of Concern (VOCs) and Variants of Interest (VOIs). This workflow has three main stages - **Genomic Analysis** , **Functional Annotation** and **Surveillance Reports** and it is being developed in combination with an interactive visualization tool [COVID-MVP](https://github.com/cidgoh/COVID-MVP). This workflow takes SARS-CoV-2 consensus sequences and Metadata (GISAID/non-GISAID) as input and produces mutation profiles which are then annotated with their respective biological functional impact using the manually curated effort [Pokay](https://github.com/nodrogluap/pokay) by Paul Gordon [@nodrogluap](https://github.com/nodrogluap).

### Genomic Analysis
This module currently supports in two different modes - "reference" & "user" which can be passed with `--mode [reference | user]`. By default, `--mode reference` is activated which allows user to build a reference library of each lineage and subsequently each variant for comparative analysis. This mode can take `fasta` file with multiple genomes (recommended & default) or single genome (`--single_genome`) with a metadata file that should have two columns (`strain`, `pango_lineage`) as minimal standard (see [Workflow Summary](#workflow-summary) for detailed options) if data is not from [GISAID](https://www.gisaid.org) (default). Data from [GISAID](https://www.gisaid.org) can be used directly by downloading from the **Genomic Epidemiology** section and passing `--gisaid` parameter.  The user mode (`--mode user`) is by default active with interactive visualization where a user can upload dataset for comparison against the reference data. Uploaded dataset can be variant called file (recommended) with `--input_type vcf` / `--input_type tsv` extensions or a `--input_type fasta` file with multiple or single genomes.

Based on the `--mode` different workflows are developed and further `--input type` determines entrance points to these workflows. The overall workflow consists of Quality Control of consensus sequences (details below in Workflow Summary), mapping consensus sequences to SARS-CoV-2 reference strain [MN908947.3 - Wuhan-Hu-1](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3), variant calling, post processing that includes flagging problematics sites; functional annotation; and mature peptide annotation. Final part of the worklfow is to collate different lineage files from the same variant and produce a surveillance report for each e.g. Delta variant report.

### Functional Annotation

In this module, the consensus sequences for each lineage is converted to a GVF (Genomic Variant Format) file and annotated with functional information.  GVF is a variant of GFF3 that is standardized for describing genomic mutations; it is used here because it can describe mutations across multiple rows, and because the "#attributes" column can store information in custom key-value pairs.  The key-value pairs added at this stage include for each mutation: VOC/VOI status, clade-defining status (for reference lineages), and functional annotations parsed from Paul Gordon's Pokay repository.  The annotated GVFs are the input to the visualization module.

### Surveillance Reports

A TSV file can be produced that summarizes mutation information for SARS-CoV-2 variants.  This file contains much of the same information found in the GVF files in a more human-readable format, with all reference lineages per variant merged into one concise TSV.


The workflow is built using [Nextflow](https://www.nextflow.io)-[DSL2](https://www.nextflow.io/docs/latest/dsl2.html), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It can use Conda/Docker/Singularity containers making installation trivial and results highly reproducible.

## Workflow Summary

The pipeline has numerous options to allow you to run only specific aspects of the workflow if you so wish. For example, for Illumina data you can skip the host read filtering step with Kraken 2 with `--skip_kraken2` or you can skip all of the assembly steps with the `--skip_assembly` parameter. See the [usage](https://nf-co.re/viralrecon/usage) and [parameter](https://nf-co.re/viralrecon/parameters) docs for all of the available options when running the pipeline.

The SRA download functionality has been removed from the pipeline (`>=2.1`) and ported to an independent workflow called [nf-core/fetchngs](https://nf-co.re/fetchngs). You can provide `--nf_core_pipeline viralrecon` when running nf-core/fetchngs to download and auto-create a samplesheet containing publicly available samples that can be accepted directly by the Illumina processing mode of nf-core/viralrecon.

### Illumina

1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Adapter trimming ([`fastp`](https://github.com/OpenGene/fastp))
4. Removal of host reads ([`Kraken 2`](http://ccb.jhu.edu/software/kraken2/); *optional*)
5. Variant calling
    1. Read alignment ([`Bowtie 2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
    2. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    3. Primer sequence removal ([`iVar`](https://github.com/andersen-lab/ivar); *amplicon data only*)
    4. Duplicate read marking ([`picard`](https://broadinstitute.github.io/picard/); *optional*)
    5. Alignment-level QC ([`picard`](https://broadinstitute.github.io/picard/), [`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    6. Genome-wide and amplicon coverage QC plots ([`mosdepth`](https://github.com/brentp/mosdepth/))
    7. Choice of multiple variant calling and consensus sequence generation routes ([`iVar variants and consensus`](https://github.com/andersen-lab/ivar); *default for amplicon data* *||* [`BCFTools`](http://samtools.github.io/bcftools/bcftools.html), [`BEDTools`](https://github.com/arq5x/bedtools2/); *default for metagenomics data*)
        * Variant annotation ([`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), [`SnpSift`](http://snpeff.sourceforge.net/SnpSift.html))
        * Consensus assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
        * Lineage analysis ([`Pangolin`](https://github.com/cov-lineages/pangolin))
        * Clade assignment, mutation calling and sequence quality checks ([`Nextclade`](https://github.com/nextstrain/nextclade))
        * Individual variant screenshots with annotation tracks ([`ASCIIGenome`](https://asciigenome.readthedocs.io/en/latest/))
    8. Intersect variants across callers ([`BCFTools`](http://samtools.github.io/bcftools/bcftools.html))
6. _De novo_ assembly
    1. Primer trimming ([`Cutadapt`](https://cutadapt.readthedocs.io/en/stable/guide.html); *amplicon data only*)
    2. Choice of multiple assembly tools ([`SPAdes`](http://cab.spbu.ru/software/spades/) *||* [`Unicycler`](https://github.com/rrwick/Unicycler) *||* [`minia`](https://github.com/GATB/minia))
        * Blast to reference genome ([`blastn`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch))
        * Contiguate assembly ([`ABACAS`](https://www.sanger.ac.uk/science/tools/pagit))
        * Assembly report ([`PlasmidID`](https://github.com/BU-ISCIII/plasmidID))
        * Assembly assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
7. Present QC and visualisation for raw read, alignment, assembly and variant calling results ([`MultiQC`](http://multiqc.info/))

### Nanopore

1. Sequencing QC ([`pycoQC`](https://github.com/a-slide/pycoQC))
2. Aggregate pre-demultiplexed reads from MinKNOW/Guppy ([`artic guppyplex`](https://artic.readthedocs.io/en/latest/commands/))
3. Read QC ([`NanoPlot`](https://github.com/wdecoster/NanoPlot))
4. Align reads, call variants and generate consensus sequence ([`artic minion`](https://artic.readthedocs.io/en/latest/commands/))
5. Remove unmapped reads and obtain alignment metrics ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
6. Genome-wide and amplicon coverage QC plots ([`mosdepth`](https://github.com/brentp/mosdepth/))
7. Downstream variant analysis:
    * Count metrics ([`BCFTools`](http://samtools.github.io/bcftools/bcftools.html))
    * Variant annotation ([`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), [`SnpSift`](http://snpeff.sourceforge.net/SnpSift.html))
    * Consensus assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
    * Lineage analysis ([`Pangolin`](https://github.com/cov-lineages/pangolin))
    * Clade assignment, mutation calling and sequence quality checks ([`Nextclade`](https://github.com/nextstrain/nextclade))
    * Individual variant screenshots with annotation tracks ([`ASCIIGenome`](https://asciigenome.readthedocs.io/en/latest/))
8. Present QC, visualisation and custom reporting for sequencing, raw reads, alignment and variant calling results ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run nf-core/viralrecon -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

4. Start running your own analysis!

    * Typical command for Illumina shotgun analysis:

        ```bash
        nextflow run nf-core/viralrecon \
            --input samplesheet.csv \
            --platform illumina \
            --protocol metagenomic \
            --genome 'MN908947.3' \
            -profile <docker/singularity/podman/conda/institute>
        ```

    * Typical command for Illumina amplicon analysis:

        ```bash
        nextflow run nf-core/viralrecon \
            --input samplesheet.csv \
            --platform illumina \
            --protocol amplicon \
            --genome 'MN908947.3' \
            --primer_set artic \
            --primer_set_version 3 \
            --skip_assembly \
            -profile <docker/singularity/podman/conda/institute>
        ```

    * Typical command for Nanopore amplicon analysis:

        ```bash
        nextflow run nf-core/viralrecon \
            --input samplesheet.csv \
            --platform nanopore \
            --genome 'MN908947.3' \
            --primer_set_version 3 \
            --fastq_dir fastq_pass/ \
            --fast5_dir fast5_pass/ \
            --sequencing_summary sequencing_summary.txt \
            -profile <docker/singularity/podman/conda/institute>
        ```

    * An executable Python script called [`fastq_dir_to_samplesheet.py`](https://github.com/nf-core/viralrecon/blob/master/bin/fastq_dir_to_samplesheet.py) has been provided if you are using `--platform illumina` and would like to auto-create an input samplesheet based on a directory containing FastQ files **before** you run the pipeline (requires Python 3 installed locally) e.g.

        ```console
        wget -L https://raw.githubusercontent.com/nf-core/viralrecon/master/bin/fastq_dir_to_samplesheet.py
        ./fastq_dir_to_samplesheet.py <FASTQ_DIR> samplesheet.csv
        ```

    * You can find the default keys used to specify `--genome` in the [genomes config file](https://github.com/nf-core/configs/blob/master/conf/pipeline/viralrecon/genomes.config). Where possible we are trying to collate links and settings for standard primer sets to make it easier to run the pipeline with standard keys; see [usage docs](https://nf-co.re/viralrecon/usage#illumina-primer-sets).

## Documentation

The nf-core/viralrecon pipeline comes with documentation about the pipeline [usage](https://nf-co.re/viralrecon/usage), [parameters](https://nf-co.re/viralrecon/parameters) and [output](https://nf-co.re/viralrecon/output).

## Acknowledgments

This workflow and scripts are written and conceptually designed by
| Name                                                      | Affiliation                                                                           |
|-----------------------------------------------------------|---------------------------------------------------------------------------------------|
| [Zohaib Anwar](https://github.com/anwarMZ)               | [Centre for Infectious Disease Genomics and One Health, Simon Fraser University, Canada](https://cidgoh.ca)                           |
| [Madeline Iseminger](https://github.com/miseminger)         | [Centre for Infectious Disease Genomics and One Health, Simon Fraser University, Canada](https://cidgoh.ca)                 |
| [Anoosha Sehar](https://github.com/Anoosha-Sehar)              | [Centre for Infectious Disease Genomics and One Health, Simon Fraser University, Canada](https://cidgoh.ca)                            |
| [Ivan Gill](https://github.com/ivansg44)                       | [Centre for Infectious Disease Genomics and One Health, Simon Fraser University, Canada](https://cidgoh.ca)                 |
| [William Hsiao](https://github.com/wwhsiao)              | [Centre for Infectious Disease Genomics and One Health, Simon Fraser University, Canada](https://cidgoh.ca)                         |
| [Paul Gordon](https://github.com/nodrogluap)                | [CSM Center for Health Genomics and Informatics, University of Calgary, Canada](http://www.ucalgary.ca/~gordonp)                         |
| [Gary Van Domselaar](https://github.com/phac-nml)                | [Public Health Agency of Canada, Canada](https://umanitoba.ca/faculties/health_sciences/medicine/units/medical_microbiology/faculty/vandomselaar.html)                         |


Many thanks to others who have helped out and contributed along the way too, including (but not limited to)\*: [Canadian COVID Genomics Network - VirusSeq, Data Analytics Working Group](https://virusseq.ca/about/governance/)


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).
For further information or help, don't hesitate to get in touch at <mzanwar@sfu.ca> or <wwhsiao@sfu.ca>

## Citations
An extensive list of references for the tools used by the workflow can be found in the [`CITATIONS.md`](CITATIONS.md) file.
