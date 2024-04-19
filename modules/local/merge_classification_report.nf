process MERGE_CLASSIFFICATION_METADATA{
    tag "$meta.id"

    conda "bioconda::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
     'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : ''}"

    input:
        tuple val(meta), path(metadata)
        tuple val(meta), path(classifier_report)

    output:
        tuple val(meta), path("Metadata_lineage.tsv"), emit: tsv

    script:
        """
        merge_classification_metadata.py --metadata ${metadata} --classifier ${pangolin_report} --output Metadata_lineage.tsv

        """
}