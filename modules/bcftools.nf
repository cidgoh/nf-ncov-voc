process BCFTOOLS{
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'
  //publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'

  tag { "BCFTOOLS_${x}" }

  input:
      tuple(path(variants), path(ref))

  output:
      path("*.vcf"), emit: vcf
  //    path("*.txt"), emit: ids

  script:
  """
  bcftools norm -f ${ref} ${variants} > ${variants.baseName}.normalized.vcf

  """
}
