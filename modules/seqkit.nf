process SEQKIT {
  //publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'

  tag { "${ids.baseName}" }

  input:
      tuple(path(ids),path(sequence))
      //path(sequence)
      //each x

  output:
      //path("*.tsv")
      path("*.fasta"), emit: fasta

  script:
    """
    seqkit grep -n -f ${ids} ${sequence} -j ${task.cpus} -o ${ids.baseName}.fasta -w 0
    """
}
