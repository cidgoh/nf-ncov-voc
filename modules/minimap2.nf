process MINIMAP2 {
    /**
    * Maps trimmed paired fastq using Minimap2 (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input
    * @output
    */

    tag {"${seq.baseName.replace("_qc", "")}"}
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.sorted.bam", mode: 'copy'

    cpus 8

    input:
        tuple(path(seq),path(ref))
        //path(indexed_reference)

    output:
        path("*.sorted.bam"), emit: bam
        path("*.bai"), emit: index

    script:
      """
      minimap2 -t ${task.cpus} -ax asm5 -a ${ref} ${seq} | \
      samtools sort -o ${seq.baseName.replace("_qc","")}.sorted.bam
      samtools index -b -@ ${task.cpus} ${seq.baseName.replace("_qc","")}.sorted.bam ${seq.baseName.replace("_qc","")}.sorted.bam.bai
      """
}
