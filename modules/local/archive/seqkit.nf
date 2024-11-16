process SEQKIT {
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":", "_")}", pattern: "*.fasta", mode: 'copy'
    tag { "${ids.Name}" }
    label 'dev_env'

    input:
    tuple val(meta), path(ids)
    tuple val(meta2), path(sequence)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta

    script:
    """
    seqkit grep -n -f ${ids} ${sequence} -j ${task.cpus} -o ${ids.baseName}.fasta -w 0
    """
}
