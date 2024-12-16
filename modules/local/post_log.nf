process POST_LOG {
  tag "${meta.id}"
  conda "conda-forge::pandas=1.4.3"
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/requests:2.26.0'
    : ''}"

  input:
  tuple val(meta), path(log_file)
  path config

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def log_date = "${params.end_date}"

  """

    virusmvp_web_connector.py \\
        -t ${params.prefix} \\
        -c ${log_file} \\
        -r ${config}
        
  """
}
