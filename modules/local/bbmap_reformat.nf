process BBMAP {

  tag { "${sequence.baseName}.fasta" }

  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'

  label 'dev_env'

  input:
      path(sequence)

  output:
      path("*.fasta"), emit: qcfasta

  when:
      sequence.size() > 0

  script:
    """
    reformat.sh \
    in=${sequence} \
    out=${sequence.baseName}.qc.fasta \
    maxns=${params.maxns} minlength=${params.minlength} addunderscore tossjunk
    """
}
