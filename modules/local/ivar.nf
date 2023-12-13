process IVAR {
    tag { "${bam.baseName}" }

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.variants.tsv", mode: 'copy'

    label 'dev_env'

    input:
    tuple(path(bam), path(ref), path(ref_gff))

    output:
    path("*.variants.tsv"), emit: variants


    script:
        """
        samtools mpileup -aa -A -d ${params.var_MaxDepth} -s --output-MQ -O --reference ${ref} ${bam} |\
        ivar variants \
        -r ${ref} \
        -m ${params.var_MinDepth} \
        -p ${bam.baseName}.variants \
        -q ${params.var_MinVariantQuality} \
        -t ${params.var_MinFreqThreshold} \
        -g ${ref_gff}
        """
}
