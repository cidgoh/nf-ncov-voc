process SNPEFF {

    tag {"snpEff${vcf.baseName}"}
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.annotated.vcf", mode: 'copy'

    cpus 1

    input:
        path(vcf)

    output:
        path("*.annotated.vcf"), emit: annotated_vcf

    script:
      """
      JAVA -jar $baseDir/snpEff/snpEff.jar \
      ann MN908947.3 -v \
      -formatEff \
      -hgvs1LetterAa \
      -hgvsOld \
      -noShiftHgvs \
      -sequenceOntology \
      -config $baseDir/snpEff/snpeff.config \
      -dataDir $baseDir/snpEff/SnpEffDB/ \
      ${vcf} \
      -htmlStats ${vcf.baseName}.snpEff.html \
      -csvStats ${vcf.baseName}.snpEff.csv \
      > ${vcf.baseName}.annotated.vcf
      """
}
