process EXTRACTVARIANTS {
      tag "$meta.id"

      conda "bioconda::pandas=1.4.3"
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : ''}"

      publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'

      input:
            tuple val(meta), path(variants)
            tuple val(meta2), path(metadata)
            val file
        
      output:
            tuple val(meta2), path("*.txt"), emit: txt

      script:

      def args = task.ext.args ?: ''
      def prefix = task.ext.prefix ?: "${meta2.id}"
      def variant_file = file ? "--variants $variants" : ''
     

      """
      parse_variants.py \\
      $variant_file \\
      --metadata ${metadata} \\
      --outfile ${prefix}_variants.txt
      """
}