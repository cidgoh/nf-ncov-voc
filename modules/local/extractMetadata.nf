process extractMetadata {
    tag "$meta.id"

    conda "bioconda::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : ''}"

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'

    input:
        tuple val(meta), path(metadata)
        tuple val(meta2), val(x)
        val time
    
    output:
    
        tuple val(meta2), path("*.tsv.gz"), emit: tsv
        tuple val(meta2), path("*.txt"), emit: txt

    script:
        def prefix = task.ext.prefix ?: "${meta2.id}"
        def time = time ? "--start_date ${params.start_date} --end_date ${params.end_date} " : ''
        """

        extract_metadata.py \\
        $time \\
        --group ${x} \\
        --table ${metadata} \\
        --outtable ${prefix}_metadata.tsv.gz \\
        --outids ${prefix}.txt
        """
}