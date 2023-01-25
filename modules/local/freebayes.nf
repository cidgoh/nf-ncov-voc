process FREEBAYES {

    tag {"${bam.baseName}"}

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'

    label 'dev_env'

    input:
    tuple(path(bam), path(ref), path(index))
    path(bam_index)

    output:
    path("*.vcf"), emit: vcf

    script:

        """
        freebayes \
        -p ${params.ploidy} \
        -f ${ref} \
        -F ${params.var_MinFreqThreshold} \
        --min-coverage ${params.var_MinDepth} \
        --pooled-continuous \
        ${bam} -v ${bam.baseName}.vcf

        """

}
