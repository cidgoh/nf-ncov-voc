process FREEBAYES {

    tag { "VariantCalling_FREEBAYES_${bam}" }
    cpus 1
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'

    input:
    tuple(path(bam), path(ref), path(index))
    path(bam_index)

    output:
    path("${bam.baseName}.vcf"), emit: variants

    script:
        """
        #freebayes-parallel \
        #          <(fasta_generate_regions.py ${index} 10000) ${task.cpus} \
        #          -p 1 \
        #          -f ${ref} \
        #          -F 0.2 \
        #          -C 1 \
        #          --pooled-continuous \
        #          --min-coverage ${params.var_MinDepth} \
        #          --gvcf --gvcf-dont-use-chunk true ${bam} |
        #          sed s/QR,Number=1,Type=Integer/QR,Number=1,Type=Float/ > ${bam.baseName}.gvcf


        freebayes \
                  -p 1 \
                  -f ${ref} \
                  -F 0.2 \
                  -C 1 \
                  --pooled-continuous \
                  --min-coverage ${params.var_MinDepth} \
                  ${bam} |
                  sed s/QR,Number=1,Type=Integer/QR,Number=1,Type=Float/ > ${bam.baseName}.vcf
        """
}
