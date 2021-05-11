process MINIMAP2 {
    /**
    * Maps trimmed paired fastq using Minimap2 (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input
    * @output
    */

    tag {"Mapping_${seq}"}
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.sorted.bam", mode: 'copy'

    cpus 4

    input:
        tuple(path(seq),path(ref))
        //path(indexed_reference)

    output:
        path("*.sorted.bam"), emit: bam
        path("*.bai"), emit: index

    script:
      """
      minimap2 -t ${task.cpus} -ax asm5 -a ${ref} ${seq} | \
      samtools sort -o ${seq.baseName}.sorted.bam
      samtools index -b -@ ${task.cpus} ${seq.baseName}.sorted.bam ${seq.baseName}.sorted.bam.bai
      """
}
