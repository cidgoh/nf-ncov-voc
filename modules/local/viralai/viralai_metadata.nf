process DOWNLOAD_VIRALAI_METADATA {
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cidgoh/virus-mvp-viralai:latest': '' }"

    
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.gz", mode: 'copy'

    output:
        path("*.gz"), emit: csv

    script:
    def prefix = "virusseq-metadata" 
    
    """
        #!/usr/bin/env bash
        dnastack config set collections.url "${params.collections_api_url}"
        dnastack collections query ${params.viralai_collection_slug_name} "SELECT * FROM collections.${params.viralai_collection_slug_name}.samples ORDER BY sample_collection_date" --format csv | gzip > ${prefix}.csv.gz
    
    """

}
