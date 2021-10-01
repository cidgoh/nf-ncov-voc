process BWAINDEX {

    tag { "indexing_SARS-CoV-2_genome" }

    input:
      path(ref)

    output:
      file("*.fasta*")

    script:
      """
      bwa index -a bwtsw $ref
      """
}
