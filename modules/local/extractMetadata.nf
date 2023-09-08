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
        def voc = x ? "--voc $x --outtable ${meta2.id}_metadata.tsv.gz --outids ${meta2.id}.txt" : ''
        //def prefix = voc ? "${meta2.id}" : "${params.enddate}"
        //def prefix = task.ext.prefix ?: "${meta2.id}"
        def time = time ? "--startdate ${params.startdate} --enddate ${params.enddate}" : ''


        """
        extract_metadata.py \\
        ${time} \\
        --table ${metadata} \\
        --criteria ${params.grouping_criteria} \\
        ${voc} 

        
        """
}