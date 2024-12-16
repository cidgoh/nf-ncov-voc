process NCOVSPLITMUTATIONSPOKAY {

  tag "${meta.id}"

  conda "conda-forge::pandas=1.4.3"
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/pandas:1.4.3'
    : 'amancevice/pandas:1.4.3'}"

  input:
  tuple val(meta), path(annotations)
  path tsv

  output:
  tuple val(meta), path("*.tsv"), emit: tsv

  script:

  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"

  """
    splitmutationnames_functionalannotation.py --functional_annotations ${annotations} \\
      --names_to_split ${tsv} \\
      ${args} \\
      --out_functions ${prefix}.processed.tsv

  """
}
