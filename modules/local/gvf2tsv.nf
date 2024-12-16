process GVF2TSV {
    tag { "${meta.id}" }
    label 'process_low'

    conda "conda-forge::pandas=1.4.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:1.4.3'
        : 'amancevice/pandas:1.4.3'}"

    input:
    tuple val(meta), path(gvf)

    output:
    tuple val(meta), path("*.tsv"), emit: surveillancetsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //def mode_command = meta.mode == 'user' ? "--user" : "--all_variants"
    """
    gvf2tsv.py \\
        --gvf_file ${gvf} \\
        ${args} \\
        --outtsv ${prefix}_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gvf2tsv: \$(gvf2tsv.py --version 2>&1 | sed 's/^.*gvf2tsv.py //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gvf2tsv: \$(gvf2tsv.py --version 2>&1 | sed 's/^.*gvf2tsv.py //; s/ .*\$//')
    END_VERSIONS
    """
}
