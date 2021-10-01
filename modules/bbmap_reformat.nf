process BBMAP {

  tag { "${sequence.baseName}.fasta" }

  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'

  input:
      path(sequence)

  output:
      path("*.fasta"), emit: qcfasta

  script:
    """
    reformat.sh \
    in=${sequence} \
    out=${sequence.baseName}.qc.fasta \
    maxns=580 minavgquality=20 minlength=29000 addunderscore tossjunk
    """
}
