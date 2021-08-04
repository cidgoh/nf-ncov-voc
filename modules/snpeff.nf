process SNPEFF {

    tag {"${filtered_vcf.baseName}"}
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.annotated.vcf", mode: 'copy'

    cpus 1

    input:
        path(filtered_vcf)

    output:
        path("*.annotated.vcf")

    script:
      """
      snpEff MN908947.3 -v \
      -formatEff \
      -hgvs1LetterAa \
      -hgvsOld \
      -noShiftHgvs \
      -sequenceOntology \
      ${filtered_vcf} \
      > ${filtered_vcf.baseName}.annotated.vcf
      """
}
