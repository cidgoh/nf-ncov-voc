process METADATA_HARMONIZER {
    tag "${meta.id}"

    conda "conda-forge::pandas=1.4.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'community.wave.seqera.io/library/pandas_pip_pyyaml:aab0deea3e114255'
        : 'community.wave.seqera.io/library/pandas_pip_pyyaml:aab0deea3e114255'}"

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
