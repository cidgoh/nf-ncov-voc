process TAGPROBLEMATICSITES_NCOV {

    tag "$meta.id"

    conda "bioconda::cyvcf=0.8.0 bioconda::gffutils=0.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cyvcf2:0.8.0--py36_0' :
        'quay.io/biocontainers/cyvcf2:0.8.0--py36_0' }"
    
    input:
        tuple val(meta), path(vcf)
        tuple val(meta2), path(prob_vcf)

    output:
        tuple val(meta), path("*.vcf"), emit: vcf

    when:
      vcf.size() > 0

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

      """
      problematic_sites_tag.py \\
        --vcffile $vcf \\
        --filter_vcf $prob_vcf \\
        --output_vcf ${prefix}.filtered.vcf \\
        $args
      """
}
