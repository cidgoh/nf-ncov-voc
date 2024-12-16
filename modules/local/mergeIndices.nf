process MERGE_INDICES {
  tag "${meta.id}"
  cpus 4

  conda "conda-forge::dask=2023.10.1"
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'biocontainers/dask:2023.10.1-py11-ol9_cv1'
    : 'daskdev/dask:2023.10.1-py3.11'}"

  input:
  tuple val(meta), path(index)
  tuple val(meta2), path(previndex)

  output:
  tuple val(meta), path("*.tsv"), emit: updated_index

  script:
  def end_date = "${params.end_date}"
  def last_index = previndex ? "--original_index ${previndex}" : ''

  """
    merge_indices.py \\
        --gvf_indices ${index} \\
        ${last_index} \\
        --cpus ${task.cpus} \\
        --index_savefile ${end_date}.index.tsv
  """
}
