process BWAINDEX {

    tag { "bwa_index_Covid_19" }

    input:
      path(ref)

    output:
      file("*.fasta*")

    script:
      """
      bwa index -a bwtsw $ref
      """
}
