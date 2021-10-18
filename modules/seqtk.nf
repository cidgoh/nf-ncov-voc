process SEQTK {
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'

  tag { "ExtractFasta_${ids}" }

  input:
      tuple(path(ids),path(sequence))
      //path(sequence)
      //each x

  output:
      //path("*.tsv")
      path("*.fasta"), emit: fasta

  script:
    """
    seqtk subseq ${sequence} ${ids} > ${ids.baseName}.fasta
    """
}
