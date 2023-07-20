process MERGE_PANGOLIN_METADATA{
    tag "$meta.id"

    conda "bioconda::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
     'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : ''}"

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'

    input:
        tuple val(meta), path(metadata)
        tuple val(meta), path(pangolin_report)

    output:
        tuple val(meta), path("Metadata_lineage.tsv"), emit: tsv

    script:
        """
        merge_pangolin_metadata.py --metadata ${metadata} --pangolin ${pangolin_report} --output Metadata_lineage.tsv

        """
}