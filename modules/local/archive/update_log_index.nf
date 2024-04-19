process UPDATE_INDEX_LOG {

  tag "$meta.id"

  conda "bioconda::pandas=1.4.3"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : '' }"
  

  input:
      tuple val(meta), path(gvfs)
      tuple val(meta), path(log_file)
      tuple val(meta), path(index)
      
  output:
      tuple val(meta), path("*.tsv"), emit: tsv
      tuple val(meta), path("*.log"), emit: log
      

  script:

  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def log_date = "${params.end_date}"

  """
    update_index_and_logfile.py \\
        --gvf_files ${gvfs} \\
        --partial_logfile ${log_file} \\
        --mutation_index ${index} \\
        --index_savefile ${log_date}.index.tsv \\
        --log_savefile ${log_date}.log 

  """

}