process SEQKIT {
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'

  tag { "${ids.Name}" }

  input:
      tuple(path(ids),path(sequence))

  output:
      path("*.fasta"), emit: fasta

  script:
    """
    seqkit grep -n -f ${ids} ${sequence} -j ${task.cpus} -o ${ids.baseName}.fasta -w 0
    """
}
