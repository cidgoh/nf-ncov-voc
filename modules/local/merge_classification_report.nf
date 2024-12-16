process MERGE_CLASSIFFICATION_METADATA {
    tag "${meta.id}"
    conda "conda-forge::pandas=1.4.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:1.4.3'
        : 'community.wave.seqera.io/library/pip_pandas:f80b46869e03f6ef'}"

    input:
    tuple val(meta), path(metadata)
    tuple val(meta2), path(classifier_report)

    output:
    tuple val(meta), path("Metadata_lineage.tsv"), emit: tsv

    script:
    """
        merge_classification_metadata.py --metadata ${metadata} --classifier ${classifier_report} --output Metadata_lineage.tsv
    """
}
