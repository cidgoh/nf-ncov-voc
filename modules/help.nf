def printHelp() {
  log.info"""
  Usage:
    nextflow run connor-lab/ncov2019-artic-nf -profile (singularity,docker,conda) ( --illumina | --nanpolish | --medaka ) --prefix [prefix] [workflow-options]

  Description:
    Variant Calling workflow for SARS-CoV-2 Variant of Concern (VOC) and Variant of Interest (VOI) consensus sequences to generate data for Visualization
    All options set via CLI can be set in conf directory

  Nextflow arguments (single DASH):
    -profile                  Allowed values: conda & singularity

  Mandatory workflow arguments (mutually exclusive):
    --prefix                  A (unique) string prefix for output directory for each run.
    --startdate               Start date (Submission date) to extract dataset (yyyy-mm-dd) (Default: "2021-01-01")
    --enddate                 Start date (Submission date) to extract dataset (yyyy-mm-dd) (Default: "")

  Optional:
    --outdir                  Output directory (Default: $baseDir/results)
    --ivar                    Run the ivar workflow (Default: false, use freebayes workflow)
    --ref                     Path to SARS-CoV-2 reference fasta file (Default: $baseDir/.github/data/refdb)
    --bwa                     Use BWA for mapping reads (Default: false, use Minimap2)
    --bwa_index               Path to BWA index files (Default: $baseDir/.github/data/refdb)
    --gff                     Path to annotation gff for variant consequence calling and typing. (Default: $baseDir/.github/data/features)
    --mpileupDepth            Mpileup depth (Default: unlimited)
    --var_FreqThreshold       Variant Calling frequency threshold for consensus variant (Default: 0.75)
    --var_MinDepth            Minimum coverage depth to call variant (ivar variants -m, freebayes -u Default: 10)
    --var_MinFreqThreshold    Minimum frequency threshold to call variant (ivar variants -t, Default: 0.25)
    --varMinVariantQuality    Minimum mapQ to call variant (ivar variants -q, Default: 20)

  """.stripIndent()
}
