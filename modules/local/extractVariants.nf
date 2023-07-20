process EXTRACTVARIANTS {
      tag "$meta.id"

      conda "bioconda::pandas=1.4.3"
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : ''}"

      publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'

      input:
            tuple val(meta), path(variants)
            tuple val(meta), path(metadata)
        
      output:
            tuple val(meta), path("*.txt"), emit: txt

      script:

      """
      parse_variants.py \
      --variants ${variants} \
      --metadata ${metadata} \
      --outfile Metadata_variants.txt
      """
}