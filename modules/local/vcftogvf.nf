process VCFTOGVF {

  tag "$meta.id"
  errorStrategy 'ignore'

  conda "bioconda::pandas=1.4.3"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : '' }"
  
  input:
      tuple val(meta), path(vcf)
      path stats
      val threshold
      tuple val(meta2), path(json)
      val lineage
      val wastewater
      val sampledesc
          
  output:
      tuple val(meta), path("*.gvf"), emit: gvf

  script:

  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def strain = lineage ? "--strain ${prefix}" : ''
  def stat     = stats ? "--size_stats ${stats}" : ''
  def wastewater = wastewater ? "--wastewater" : ''
  def group =  "--sample_group ${meta.id}"


  """
    vcf2gvf.py --vcffile $vcf \\
      $stat \\
      --clades_threshold $threshold \\
      --gene_positions $json \\
      $wastewater \\
      $strain \\
      --sample_desc $sampledesc \\
      $group \\
      $args \\
      --outgvf ${prefix}.gvf

  """

}