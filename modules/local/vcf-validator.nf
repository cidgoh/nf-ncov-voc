process VCF_VALIDATOR {
    tag "${meta.id}"
    conda "bioconda::vcf-validator=0.10.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/vcf-validator:0.10.0--h9cfbc0b_2'
        : 'quay.io/biocontainers/vcf-validator:0.10.0--h9cfbc0b_2'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.txt"), emit: validation_summary

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    vcf_validator -l error -r summary < ${vcf}
    """

    stub:
    """
    echo "VCF validator stub output" > summary
    mkdir -p vcf-validator-test
    touch vcf-validator-test/dummy.txt
    """
}
