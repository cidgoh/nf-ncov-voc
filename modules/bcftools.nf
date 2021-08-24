process BCFTOOLS{
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'
  tag {"${variants.baseName.replace("variants", "")}"}

  input:
      tuple(path(variants), path(ref))

  output:
      path("*.vcf"), emit: normalized_vcf

  script:
  """
  bcftools norm -f ${ref} ${variants} > ${variants.baseName}.normalized.vcf

  """
}
