process BWAMEM {

    tag {"${seq.baseName}"}

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.sorted.bam", mode: 'copy'

    label 'dev_env'

    input:
        tuple(path(seq),path(ref))
        path(indexed_reference)

    output:
        path("*.sorted.bam"), emit: bam

    when:
      seq.size() > 0

    script:
      """
      bwa mem -t ${task.cpus} ${ref} ${seq} | \
      samtools sort -o ${seq.baseName}.sorted.bam
      """
}
