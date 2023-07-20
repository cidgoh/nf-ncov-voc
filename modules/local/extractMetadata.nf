process extractMetadata {
    tag "$meta.id"

    conda "bioconda::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : ''}"

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'

    input:
        tuple val(meta), path(metadata)
        tuple val(meta), val(x)
    
    output:
    
        tuple val(meta), path("*.tsv.gz"), emit: tsv
        path("*.txt"), emit: txt
    
    when:
        x.size() > 0

    script:
        def voc = x ? "--voc $x" : ''
        """
        extract_metadata.py \\
        --startdate ${params.startdate} \\
        --enddate ${params.enddate} \\
        --table ${metadata} \\
        --criteria ${params.grouping_criteria} \\
        ${voc}
        
        """
}