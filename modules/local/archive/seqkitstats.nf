process SEQKITSTATS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.2.0--h9ee0642_0':
        'biocontainers/seqkit:2.2.0--h9ee0642_0' }"
  
    input:
        path(samples)

    output:
        path("samples_stats.tsv"), emit: stats

    script:
        """

        seqkit stats ${samples} -o samples_stats.tsv
        """
}
