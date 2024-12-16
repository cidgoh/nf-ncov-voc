process FUNCTIONALANNOTATION {
  tag "${meta.id}"
  conda "conda-forge::pandas=1.4.3"
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/pandas:1.4.3'
    : 'amancevice/pandas:1.4.3'}"

  input:
  tuple val(meta), path(gvf)
  tuple val(meta2), path(tsv)

  output:
  tuple val(meta), path("*.gvf"), emit: gvf

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"

  """
    addfunctions2gvf.py \\
      --ingvf ${gvf} \\
      --outgvf ${prefix}.annotated.gvf \\
      --functional_annotations ${tsv}

  """
}
