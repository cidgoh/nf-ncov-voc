process SNPEFF {

    tag {"${filtered_vcf.baseName}"}

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'

    label 'dev_env'

    input:
        path(filtered_vcf)

    output:
        path("*.vcf"), emit: peptide_vcf

    script:
      """
      snpEff NC_045512.2 -v \
      -formatEff \
      -hgvs1LetterAa \
      -hgvs \
      -hgvsOld \
      ${filtered_vcf} \
      > ${filtered_vcf.baseName}.SNPEFF.vcf
      """
}
