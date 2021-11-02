process FREEBAYES {

    tag {"${bam.baseName}"}

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.gvcf", mode: 'copy'

    cpus 1

    input:
    tuple(path(bam), path(ref), path(index))
    path(bam_index)

    output:
    path("*.gvcf"), emit: gvcf

    script:

        """
        freebayes \
        -p 1 \
        -f ${ref} \
        -F 0.2 \
        -C 1 \
        --pooled-continuous \
        --min-coverage ${params.var_MinDepth} \
        --limit-coverage ${params.var_downsample} \
        ${bam} |
        sed s/QR,Number=1,Type=Integer/QR,Number=1,Type=Float/ > ${bam.baseName}.gvcf

        """

}
