process EXTRACTVARIANTS {
      tag "${meta.id}"
      conda "conda-forge::pandas=1.4.3"
      container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
            ? 'https://depot.galaxyproject.org/singularity/pandas:1.4.3'
            : 'amancevice/pandas:1.4.3'}"

      input:
      tuple val(meta), path(variants)
      tuple val(meta2), path(metadata)
      val variant_file
      val virusseq
      val criteria
      val variable
      val time

      output:
      tuple val(meta), path("*.txt"), emit: txt
      tuple val(meta), path("*.log"), emit: log, optional: true

      script:
      def prefix = task.ext.prefix ?: "${meta.id}"
      def val_variant_file = variant_file ? "--variants ${variants}" : ''
      def virusseq_update = virusseq ? "--virusseq" : ''
      def val_time = time ? "--start_date ${params.start_date} --end_date ${params.end_date}" : ''
      def val_variable = variable ? "--variable ${params.variable}" : ''
      def val_criteria = criteria ? "--criteria ${criteria}" : ''

      """
      parse_variants.py \\
      ${val_variant_file} \\
      ${virusseq_update} \\
      ${val_variable} \\
      ${val_time} \\
      ${val_criteria} \\
      --metadata ${metadata} \\
      --outfile ${prefix}_variants.txt \\
      --logfile ${prefix}_web.log
      """
}
