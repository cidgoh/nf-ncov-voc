process virusseqMapLineage {
  tag { "Mapping VirusSeq dataset to GISAID for lineages" }
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":", "_")}", pattern: "*.tsv", mode: 'copy'

  input:
  tuple val(meta), path(metadata)
  tuple val(meta2), path(gisaid_metadata)

  output:
  path ("VirusSeq_Mapped.tsv"), emit: mapped

  script:
  """
      map_virusseq_GISAID.py --virusseq ${metadata} --gisaid ${gisaid_metadata} --output VirusSeq_Mapped.tsv

      """
}

process grabIndex {
  tag { "grabing_SARS-CoV-2_index" }

  input:
  path index_folder

  output:
  file "*.fasta.*"

  script:
  """
      ln -sf ${index_folder}/*.fasta.* ./
      """
}

process IVAR_VARIANTS_TO_VCF {
  tag "${meta.id}"
  conda "conda-forge::python=3.9.5 conda-forge::matplotlib=3.5.1 conda-forge::pandas=1.3.5 conda-forge::r-sys=3.4 conda-forge::regex=2021.11.10 conda-forge::scipy=1.7.3 conda-forge::biopython=1.79"
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0'
    : 'quay.io/biocontainers/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0'}"

  input:
  tuple val(meta), path(tsv)

  output:
  tuple val(meta), path("*.vcf"), emit: vcf
  tuple val(meta), path("*.log"), emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    ivar_variants_to_vcf.py \\
        ${tsv} \\
        ${prefix}.vcf \\
        ${args} \\
        > ${prefix}.variant_counts.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

process processGVCF {
  tag "${meta.id}"
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":", "_")}", pattern: "*variants.vcf", mode: 'copy'

  input:
  tuple val(meta), path(gvcf)

  output:
  tuple val(meta), path("*.vcf"), emit: vcf
  tuple val(meta), path("*.log"), emit: log
  path "versions.yml", emit: versions

  when:
  gvcf.size() > 0

  script:
  """
      process_gvcf.py -d ${params.var_MinDepth} \
      -l ${params.lower_ambiguityFrequency} \
      -u ${params.upper_ambiguityFrequency} \
      -m ${gvcf.baseName}.mask.txt \
      -v ${gvcf.baseName}.variants.vcf \
      -c ${gvcf.baseName}.consensus.vcf ${gvcf}
      """
}


process surveillanceRawTsv {
  tag { "Generating Raw Surveillance Data (TSV)" }
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":", "_")}", pattern: "*.tsv", mode: 'copy'

  input:
  path gvf
  tuple path(variants), path(stats)

  output:
  path ('*.tsv'), emit: surveillancetsv

  script:
  if (params.mode == 'user') {
    """
    gvf2tsv.py --gvf_files ${gvf} \
    --user
    """
  }
  else {
    """
    gvf2tsv.py --gvf_files ${gvf} \
    --clades ${variants} \
    --table ${stats} \
    --all_variants

    """
  }
}

process surveillancePDF {
  tag { "Generating Surveillance reports PDF" }
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":", "_")}", pattern: "*.pdf", mode: 'copy'

  input:
  each tsv
  path surveillanceindicators
  path metadata

  output:
  path ("*.pdf"), emit: surveillance_pdf

  script:
  if (params.mode == 'reference') {

    """
    surveillance_report_pdf.py --tsv ${tsv} \
    --functions_table ${surveillanceindicators} \
    --metadata ${metadata} \
    --virusseq True > ${tsv.baseName}.tex

    tectonic -X compile ${tsv.baseName}.tex --reruns 3 --keep-intermediates --keep-logs


    """
  }
  else {
    """
    surveillance_report_pdf.py --tsv ${tsv} \
    --functions_table ${surveillanceindicators} > ${tsv.baseName}.tex 
    
    tectonic -X compile ${tsv.baseName}.tex --reruns 3 --keep-intermediates --keep-logs

    """
  }
}
