process DOWNLOAD_VIRALAI_MULTIFASTA {
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://cidgoh/virus-mvp-viralai:latest'
        : 'cidgoh/virus-mvp-viralai:latest'}"

    output:
    path ("*.xz"), emit: xz

    script:
    """
        #!/usr/bin/env bash
        dnastack config set collections.url "${params.collections_api_url}"
        dnastack config set drs.url "${params.collections_drs_url}"
        dnastack collections query ${params.collection_slug_name} "SELECT * FROM viralai.${params.collection_slug_name}.files WHERE NAME LIKE 'CanCOGeN/multifasta_compressed/%' ORDER BY created_time DESC LIMIT 1" | jq -r '.[].drs_url' > test.txt
        dnastack files download -i test.txt
    """
}
