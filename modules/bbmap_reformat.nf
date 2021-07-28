process BBMAP {
  //publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.log", mode: 'copy'

  tag { "${sequence.baseName}" }

  input:
      path(sequence)

  output:
      path("*.fasta"), emit: qcfasta

  script:
    """
    reformat.sh \
    in=${sequence} \
    out=${sequence.baseName}_qc.fasta \
    maxns=145 addunderscore tossjunk > ${sequence.baseName}.log
    """
}
