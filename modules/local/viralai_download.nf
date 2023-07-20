#!/usr/bin/env nextflow

params.collection_name="virusseq"
params.collections_api_url="https://viral.ai/api/"
params.limit = "10"

process downloadFirstTenFiles {
    container 'gcr.io/dnastack-pub-container-store/dnastack-cli:latest'
    containerOptions = '--user root'

    output:
    file 'outputs/*' into downloadOutput

    script:
    """
    #!/usr/bin/env bash
    mkdir outputs
    dnastack config set collections.url "${params.collections_api_url}"
    collection_slug_name=`dnastack collections list | jq -r '.[] | select(.name == "${params.collection_name}") | .slugName'`
    query="SELECT drs_url FROM \"viralai\".\"\$collection_slug_name\".\"files\" WHERE name LIKE '%multifasta_compressed%' LIMIT ${params.limit}"
    dnastack collections query \$collection_slug_name "\$query" | jq -r '.[].drs_url' | dnastack files download -o outputs
    """
}


downloadOutput.flatMap().subscribe { println "File: ${it.name}" }