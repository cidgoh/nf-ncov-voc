# nf-ncov-voc

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)


## Introduction

**nf-ncov-voc** is a bioinformatics analysis workflow used for performing variant calling on SARS-CoV-2 genomes to identify and profile mutations in Variants of Concern (VOCs) and Variants of Interest (VOIs). This workflow has three main stages - **Genomic Analysis** , **Functional Annotation** and **Surveillance Reports**. This workflow is developed in combination with an interactive visualization tool [COVID-MVP](https://github.com/cidgoh/COVID-MVP). As an input, **nf-ncov-voc** takes SARS-CoV-2 consensus sequences and Metadata (GISAID/non-GISAID) and produces mutation profiles which are then annotated with their respective biological functional impact using the manually curated effort [Pokay](https://github.com/nodrogluap/pokay) lead by Paul Gordon [@nodrogluap](https://github.com/nodrogluap).

The workflow is built using [Nextflow](https://www.nextflow.io)-[DSL2](https://www.nextflow.io/docs/latest/dsl2.html), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It can use **Conda**/**Docker**/**Singularity** containers making installation trivial and results highly reproducible.

### Genomic Analysis

This module currently supports two different modes - "_reference_" & "_user_" which can be passed with `--mode reference` or `--mode user`. By default, `--mode reference` is activated which allows user to build a reference library for each lineage and subsequently each variant for comparative analysis. This mode can take `FASTA` file with multiple genomes (**recommended** & **default**) or single genome (passed as `--single_genome`) with a metadata file that should have two columns atleast (`strain`, `pango_lineage`) as minimal metadata (see [Workflow Summary](#workflow-summary) for detailed options). Data can directly be used from [GISAID](https://www.gisaid.org) after downloading from the **Genomic Epidemiology** section and passing `--gisaid` parameter. Similarly, _non-GISAID_ (**default**) data can also be used with minimal metadata file. The user mode (`--mode user`) is by default active when using interactive visualization through [COVID-MVP](https://github.com/cidgoh/COVID-MVP) where a user can upload dataset for comparative analysis against the reference data. Uploaded dataset can be a variant called file `VCF` or `TSV` with `--input_type vcf` (**recommended**) or `--input_type tsv`. Alternatively, a `FASTA` file can be provided using `--input_type fasta` file with a single or multiple genomes.

Based on the `--mode` different workflows are developed and further `--input type` determines entrance points to these workflows. The overall workflow consists of Quality Control of consensus sequences (details below in Workflow Summary), mapping consensus sequences to SARS-CoV-2 reference strain [MN908947.3 - Wuhan-Hu-1](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3), variant calling, post processing that includes flagging problematics sites; functional annotation; and mature peptide annotation (see [Usage](#usage) for guideline). Final part of the workflow is to collate different lineage files from the same variant and produce a surveillance report for each e.g. Delta variant report.

### Functional Annotation

In this module, the variant called `VCF` file for each lineage is converted into a `GVF` (Genomic Variant Format) file and annotated with functional information using [Pokay](https://github.com/nodrogluap/pokay). GVF is a variant of GFF3 format that is standardized for describing genomic mutations; it is used here because it can describe mutations across multiple rows, and because the "#attributes" column can store information in custom key-value pairs. The key-value pairs added at this stage include for each mutation: VOC/VOI status, clade-defining status (for reference lineages), and functional annotations parsed using [vcf2gvf.py](https://github.com/cidgoh/nf-ncov-voc/blob/master/bin/vcf2gvf.py) file written in python.

### Surveillance Reports

Different `GVF` files for the same variant are then collated and summarized into a `TSV` file that contains mutation prevalence, profile and functional impact. This feature can be used for identify and track transmission trends in a dataset, aid detection of new cluster important mutations with severe impact based on the datasets used.

### nf-ncov-voc Dataflow

![DataFlow](figs/COVIDMVP.drawio.png)

## Workflow Summary

The workflow has numerous options to allow you to run workflow with modes and alternate options for major step if you so wish. For example, in `mode --reference` user can use `BWAMEM` using `--bwa` instead of `MINIMAP2` (*default*) for mapping consensus sequences to reference genome. Similarly, `ivar` with parameter `--ivar` for variant calling instead of `freebayes` (*default*) option.

See the [parameters]() docs for all of the available options when running the workflow.

### Reference Mode

* _Data Extraction & Quality Control_
    1.  Metadata extraction ([`bin/extract_metadata.py`](https://github.com/cidgoh/nf-ncov-voc/blob/master/bin/extract_metadata.py) && [`modules/custom.nf/extractMetadata`](https://github.com/cidgoh/nf-ncov-voc/blob/master/modules/custom.nf))
    2.  Sequence extraction ([`SEQKIT`](https://github.com/shenwei356/seqkit))
    3.  Consensus QC ([`BBMAP`](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/))
* _Variant Calling_
    1.  Mapping ([`Minimap2`](); *default* [`BWA`](); *optional* )
    2.  Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    3.  Choice of multiple variant calling routes
      1.  [`Freebayes`]() *default & recommended*;
      2.  [`iVar variants`](https://github.com/andersen-lab/ivar); *optional*)
* _Post-Processing_
    1.  Filtering Problematic sites ([`problematic_sites_tag.py`](https://github.com/cidgoh/nf-ncov-voc/blob/master/bin/problematic_sites_tag.py) using [ProblematicSites_SARS-CoV-2](https://github.com/W-L/ProblematicSites_SARS-CoV2))
    2.  Variant annotation ([`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html))
    3.  Peptide annotation ([]())

### User Mode

* _Data Extraction & Quality Control_
    1.  Metadata extraction ([`bin/extract_metadata.py`](https://github.com/cidgoh/nf-ncov-voc/blob/master/bin/extract_metadata.py) && [`modules/custom.nf/extractMetadata`](https://github.com/cidgoh/nf-ncov-voc/blob/master/modules/custom.nf))
    2.  Sequence extraction ([`SEQKIT`](https://github.com/shenwei356/seqkit))
    3.  Consensus QC ([`BBMAP`](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/))
* _Variant Calling_
    1.  Mapping ([`Minimap2`](); *default* [`BWA`](); *optional* )
    2.  Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    3.  Choice of multiple variant calling routes
      1.  [`Freebayes`]() *default & recommended*;
      2.  [`iVar variants`](https://github.com/andersen-lab/ivar); *optional*)
* _Post-Processing_
    1.  Filtering Problematic sites ([`problematic_sites_tag.py`](https://github.com/cidgoh/nf-ncov-voc/blob/master/bin/problematic_sites_tag.py) using [ProblematicSites_SARS-CoV-2](https://github.com/W-L/ProblematicSites_SARS-CoV2))
    2.  Variant annotation ([`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html))
    3.  Peptide annotation ([]())


## Usage

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html) for full pipeline reproducibility _see [recipes](https://github.com/cidgoh/nf-ncov-voc/tree/master/environments)_

3. Download the pipeline and run with help for detailed parameter options:

    ```console
    nextflow run nf-ncov-voc/main.nf --help
    ```

    ```bash
    N E X T F L O W  ~  version 21.04.3
    Launching `main.nf` [berserk_austin] - revision: 93ccc86071

    Usage:
      nextflow run main.nf -profile [singularity | docker | conda) --prefix [prefix] --mode [reference | user]  [workflow-options]

    Description:
      Variant Calling workflow for SARS-CoV-2 Variant of Concern (VOC) and Variant of Interest (VOI) consensus sequences to generate data for Visualization
      All options set via CLI can be set in conf directory

    Nextflow arguments (single DASH):
      -profile                  Allowed values: conda & singularity

    Mandatory workflow arguments (mutually exclusive):
      --prefix                  A (unique) string prefix for output directory for each run.
      --mode                    A flag for user uploaded data through visualization app or high-throughput analyses (reference | user) (Default: reference)

    Optional:
      --variants                Provide a variants file (tsv) (Default: /Users/au572806/GitHub/nf-ncov-voc/.github/data/variants/variants_who.tsv)
      --input_type              Specify type of input file (vcf | tsv | fasta) (Default: vcf)
      --gisaid                  Specify if the dataset is from GISAID (gisaid) (Default: None)
      --single-genome           Specify if the dataset is single genome (single-genome) (Default: None)
      --userfile                Specify userfile (fasta | tsv | vcf) (Default: None)
      --outdir                  Output directory (Default: /Users/au572806/GitHub/nf-ncov-voc/results)
      --ivar                    Run the ivar workflow (Default: false, use freebayes workflow)
      --startdate               Start date (Submission date) to extract dataset (yyyy-mm-dd) (Default: "None")
      --enddate                 Start date (Submission date) to extract dataset (yyyy-mm-dd) (Default: "None")
      --ref                     Path to SARS-CoV-2 reference fasta file (Default: /Users/au572806/GitHub/nf-ncov-voc/.github/data/refdb)
      --bwa                     Use BWA for mapping reads (Default: false, use Minimap2)
      --bwa_index               Path to BWA index files (Default: /Users/au572806/GitHub/nf-ncov-voc/.github/data/refdb)
      --gff                     Path to annotation gff for variant consequence calling and typing. (Default: /Users/au572806/GitHub/nf-ncov-voc/.github/data/features)
      --mpileupDepth            Mpileup depth (Default: unlimited)
      --var_FreqThreshold       Variant Calling frequency threshold for consensus variant (Default: 0.75)
      --var_MinDepth            Minimum coverage depth to call variant (ivar variants -m, freebayes -u Default: 10)
      --var_MinFreqThreshold    Minimum frequency threshold to call variant (ivar variants -t, Default: 0.25)
      --varMinVariantQuality    Minimum mapQ to call variant (ivar variants -q, Default: 20)
    ```

4. Start running your own analysis!

    * Typical command for Illumina shotgun analysis:

        ```bash
        nextflow run nf-ncov-voc/main.nf \
            --input samplesheet.csv \
            --platform illumina \
            --protocol metagenomic \
            --genome 'MN908947.3' \
            -profile <docker/singularity/podman/conda/institute>
        ```

    * Typical command for Illumina amplicon analysis:

        ```bash
        nextflow run nf-ncov-voc/main.nf \
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
        nextflow run nf-ncov-voc/main.nf \
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

## Acknowledgments

This workflow and scripts are written and conceptually designed by
| Name                                                      | Affiliation                                                                           |
|-----------------------------------------------------------|---------------------------------------------------------------------------------------|
| Zohaib Anwar; [@anwarMZ](https://github.com/anwarMZ)               | [Centre for Infectious Disease Genomics and One Health, Simon Fraser University, Canada](https://cidgoh.ca)                           |
| Madeline Iseminger; [@miseminger](https://github.com/miseminger)         | [Centre for Infectious Disease Genomics and One Health, Simon Fraser University, Canada](https://cidgoh.ca)                 |
| Anoosha Sehar; [@Anoosha-Sehar](https://github.com/Anoosha-Sehar)              | [Centre for Infectious Disease Genomics and One Health, Simon Fraser University, Canada](https://cidgoh.ca)                            |
| Ivan Gill; [@ivansg44](https://github.com/ivansg44)                       | [Centre for Infectious Disease Genomics and One Health, Simon Fraser University, Canada](https://cidgoh.ca)                 |
| William Hsiao; [@wwhsiao](https://github.com/wwhsiao)              | [Centre for Infectious Disease Genomics and One Health, Simon Fraser University, Canada](https://cidgoh.ca)                         |
| Paul Gordon; [@nodrogluap](https://github.com/nodrogluap)                | [CSM Center for Health Genomics and Informatics, University of Calgary, Canada](http://www.ucalgary.ca/~gordonp)                         |
| Gary Van Domselaar; [@phac-nml](https://github.com/phac-nml)                | [Public Health Agency of Canada](https://umanitoba.ca/faculties/health_sciences/medicine/units/medical_microbiology/faculty/vandomselaar.html)                         |

Many thanks to others who have helped out and contributed along the way too, including (but not limited to)\*: [Canadian COVID Genomics Network - VirusSeq, Data Analytics Working Group](https://virusseq.ca/about/governance/)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).
For further information or help, don't hesitate to get in touch at <mzanwar@sfu.ca> or <wwhsiao@sfu.ca>

## Citations
An extensive list of references for the tools used by the workflow can be found in the [`CITATIONS.md`](CITATIONS.md) file.
