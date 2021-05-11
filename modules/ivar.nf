process IVAR {

    tag { "VariantCalling_IVAR_${bam}" }

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.variants.tsv", mode: 'copy'

    input:
    tuple(path(bam), path(ref), path(ref_gff))

    output:
    tuple path("${bam.baseName}.variants.tsv")


    script:
        """
        samtools mpileup -aa -A -d 0 --reference ${ref} -B -Q 13 ${bam} |\
        ivar variants \
        -r ${ref} \
        -m ${params.ivarMinDepth} \
        -p ${bam.baseName}.variants \
        -q ${params.ivarMinVariantQuality} \
        -t ${params.ivarMinFreqThreshold} \
        -g ${ref_gff}
        """
}
