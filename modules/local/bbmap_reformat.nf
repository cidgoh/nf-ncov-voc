process BBMAP {
    tag "$meta.id"

    conda "bioconda::bbmap=30.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0' : ''}"
    

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'

    label 'dev_env'

    input:
        tuple val(meta), path(sequence)

    output:
        tuple val(meta), path("*.fasta"), emit: fasta

    when:
        sequence.size() > 0

    script:
    def args = task.ext.args ?: ''
    """
    reformat.sh \\
    in=${sequence} \\
    out=${sequence.baseName}.qc.fasta \\
    
    """
}
