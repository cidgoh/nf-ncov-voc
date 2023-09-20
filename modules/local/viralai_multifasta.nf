process DOWNLOAD_VIRALAI_MULTIFASTA {
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cidgoh/virus-mvp-viralai:latest': '' }"

    
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.xz", mode: 'copy'

    output:
        path("*.xz"), emit: xz
    
    """
        #!/usr/bin/env bash
        dnastack config set collections.url "${params.collections_api_url}"
        dnastack config set drs.url "${params.collections_drs_url}"
        dnastack collections query ${params.viralai_collection_slug_name} "SELECT * FROM viralai.${params.viralai_collection_slug_name}.files WHERE NAME LIKE '%multifasta_compressed%' ORDER BY created_time DESC LIMIT 1" | jq -r '.[].drs_url' | dnastack files download
    
    """

}
