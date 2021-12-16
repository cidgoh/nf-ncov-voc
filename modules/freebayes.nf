process FREEBAYES {

    tag {"${bam.baseName}"}

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.gvcf", mode: 'copy'

    label 'dev_env'

    input:
    tuple(path(bam), path(ref), path(index))
    path(bam_index)

    output:
    path("*.gvcf"), emit: gvcf

    script:

        """
        freebayes \
        -p ${params.ploidy} \
        -f ${ref} \
        -F ${params.var_MinFreqThreshold} \
        -C ${params.var_MinDepth} \
        --pooled-continuous \
        --min-coverage ${params.var_MinDepth} \
        ${bam} |
        sed s/QR,Number=1,Type=Integer/QR,Number=1,Type=Float/ > ${bam.baseName}.gvcf

        """

}
