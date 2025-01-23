process TSV2PDF {
    tag "${meta.id}"
    conda "conda-forge::pandas=1.4.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'community.wave.seqera.io/library/pip_pandas_reportlab:14a90fc25f86e9ba'
        : 'community.wave.seqera.io/library/pip_pandas_reportlab:14a90fc25f86e9ba'}"

    input:
    tuple val(meta), path(tsv), path(metadata)
    tuple val(meta2), path(surveillanceindicators)

    output:
    tuple val(meta), path("*.pdf"), emit: surveillance_pdf

    script:
    def metadata_arg = metadata ? "--metadata ${metadata}" : ""
    def virusseq_arg = params.mode == 'reference' ? "--virusseq" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    surveillance_report_pdf.py \
        --tsv ${tsv} \
        --functions_table ${surveillanceindicators} \
        ${metadata_arg} \
        ${virusseq_arg} \
        --output ${prefix}_annotated.pdf
    """
}
