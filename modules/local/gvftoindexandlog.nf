process GVF_TO_INDEX_LOG {
  tag "${meta.id}"
  conda "conda-forge::dask=2023.10.1"
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'biocontainers/dask:2023.10.1-py11-ol9_cv1'
    : 'quay.io/biocontainers/dask:2023.10.1-py11-ol9_cv1'}"

  input:
  tuple val(meta), path(gvf)

  output:
  tuple val(meta), path("*.tsv"), emit: index
  tuple val(meta), path("*.log"), emit: log

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def end_date = "${params.end_date}"

  """
    gvf2indexandlog.py \\
        --gvf_file ${gvf} \\
        --index_savefile ${prefix}.${end_date}.index.tsv \\
        --log_savefile ${prefix}.${end_date}.log
  """
}
