process METADATA_HARMONIZER {
    tag "${meta.id}"

    conda "conda-forge::pandas=1.4.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:1.4.3'
        : 'amancevice/pandas:1.4.3'}"

    input:
    tuple val(meta), path(metadata)
    path yaml
    val dataset

    output:
    tuple val(meta), path("*.csv.gz"), emit: gz

    script:

    """
        harmonize_metadata.py \\
        --metadata ${metadata} \\
        --config ${yaml} \\
        --data ${dataset} \\
        --outfile harmonized_metadata.csv.gz

        """
}
