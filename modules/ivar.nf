process IVAR {

    tag { "VariantCalling_IVAR_${bam}" }

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.variants.tsv", mode: 'copy'

    input:
    tuple(path(bam), path(ref), path(ref_gff))

    output:
    tuple path("${bam.baseName}.variants.tsv")


    script:
        """
        samtools mpileup -aa -A -d 0 --reference ${ref} -B -Q 0 ${bam} |\
        ivar variants \
        -r ${ref} \
        -m ${params.var_MinDepth} \
        -p ${bam.baseName}.variants \
        -q ${params.var_MinVariantQuality} \
        -t ${params.var_MinFreqThreshold} \
        -g ${ref_gff}
        """
}
