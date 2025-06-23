process PATHOPLEXUS_METADATA {
    tag "Fetch MPOX metadata"
    label 'process_low'

    conda "conda-forge::pandas=2.2.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.3'
        : 'quay.io/biocontainers/pandas:2.2.3'}"

    input:
    val output_filename
    val limit

    output:
    path "${output_filename}", emit: metadata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    pathoplexus_metadata.py \\
        -o ${output_filename} \\
        -l ${limit} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${output_filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
