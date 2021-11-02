process BBMAP {

  tag { "${sequence.baseName}.fasta" }

  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'

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
    maxns=580 minlength=29000 addunderscore tossjunk
    """
}
