process SEQKITSTATS {
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'

  tag { "Data Statistics" }

  input:
      path(samples)

  output:
      path("samples_stats.tsv"), emit: stats

  script:
    """
    seqkit stats ${samples} -o samples_stats.tsv
    """
}
