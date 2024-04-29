process PROCESS_VIRALAI_METADATA {
    tag "$meta.id"

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : ''}"

    input:
        tuple val(meta), path(metadata)
        path(alias)
        
    output:
        tuple val(meta), path("*.csv.gz"), emit: gz


    script:

        """
        process_viralaimetadata.py \\
        --csv $metadata \\
        --alias $alias \\
        --outfile viralai_metadata_processed.csv.gz

        """
}    
