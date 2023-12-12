process EXTRACTVARIANTS {
      tag "$meta.id"

      conda "bioconda::pandas=1.4.3"
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : ''}"

      publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.txt", mode: 'copy'

      input:
            tuple val(meta), path(variants)
            tuple val(meta2), path(metadata)
            val variant_file
            val virusseq
            val criteria
            val variable
            val time
        
      output:
            tuple val(meta2), path("*.txt"), emit: txt
            tuple val(meta2), path("*.log"), emit: log, optional: true

      script:

      def args = task.ext.args ?: ''
      def prefix = task.ext.prefix ?: "${meta2.id}"
      def variant_file = variant_file ? "--variants $variants" : ''
      def virusseq_update = virusseq ? "--virusseq" : ''
      def time = time ? "--start_date ${params.start_date} --end_date ${params.end_date}" : ''
      def variable = variable ? "--variable ${params.variable}" : ''
      def criteria = criteria ? "--criteria $criteria" : ''

      """
      parse_variants.py \\
      $variant_file \\
      $virusseq_update \\
      $variable \\
      $time \\
      $criteria \\
      --metadata ${metadata} \\
      --outfile ${prefix}_variants.txt \\
      --logfile ${prefix}_web.log
      """
}