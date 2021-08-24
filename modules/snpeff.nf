process SNPEFF {

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'
    tag {"${filtered_vcf.baseName.replace(".variants.normalized.filtered", "")}"}

    input:
        path(filtered_vcf)

    output:
        path("*.vcf"), emit: peptide_vcf

    script:
      """
      snpEff MN908947.3 -v \
      -formatEff \
      -hgvs1LetterAa \
      -hgvsOld \
      -noShiftHgvs \
      -sequenceOntology \
      ${filtered_vcf} \
      > ${filtered_vcf.baseName}.SNPEFF_annotated.vcf
      """
}
